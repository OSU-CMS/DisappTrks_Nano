"""Numerical core of the electron and muon fiducial-map construction."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np


@dataclass(frozen=True)
class HotSpot:
    eta: float
    phi: float
    radius: float
    sigma: float


@dataclass(frozen=True)
class FiducialMapSummary:
    mean_inefficiency: float
    stddev_inefficiency: float
    inefficiency: np.ndarray
    significance: np.ndarray
    hot_spots: tuple[HotSpot, ...]

    def as_dict(self) -> dict[str, Any]:
        return {
            "mean_inefficiency": self.mean_inefficiency,
            "stddev_inefficiency": self.stddev_inefficiency,
            "hot_spots": [asdict(hot_spot) for hot_spot in self.hot_spots],
        }


def summarize_fiducial_map(
    before: np.ndarray,
    after: np.ndarray,
    eta_edges: np.ndarray,
    phi_edges: np.ndarray,
    threshold: float = 2.0,
) -> FiducialMapSummary:
    """Reproduce the legacy bin-by-bin inefficiency hot-spot calculation."""
    before = np.asarray(before, dtype=float)
    after = np.asarray(after, dtype=float)
    eta_edges = np.asarray(eta_edges, dtype=float)
    phi_edges = np.asarray(phi_edges, dtype=float)

    if before.shape != after.shape:
        raise ValueError("before and after histograms must have identical shapes")
    if before.shape != (len(eta_edges) - 1, len(phi_edges) - 1):
        raise ValueError("histogram shape does not match the supplied axes")

    occupied = before > 0.0
    if not np.any(occupied):
        raise ValueError("before histogram has no occupied bins")

    mean = float(after[occupied].sum() / before[occupied].sum())
    inefficiency = np.zeros_like(after)
    inefficiency[occupied] = after[occupied] / before[occupied]

    n_occupied = int(np.count_nonzero(occupied))
    stddev = (
        float(
            np.sqrt(
                np.sum((inefficiency[occupied] - mean) ** 2)
                / (n_occupied - 1)
            )
        )
        if n_occupied > 1
        else 0.0
    )
    significance = np.zeros_like(after)
    if stddev > 0.0:
        significance[occupied] = (inefficiency[occupied] - mean) / stddev

    hot_spots = []
    for ix, iy in np.argwhere(
        occupied & (inefficiency > 0.0) & (significance > threshold)
    ):
        eta_width = eta_edges[ix + 1] - eta_edges[ix]
        phi_width = phi_edges[iy + 1] - phi_edges[iy]
        hot_spots.append(
            HotSpot(
                eta=float((eta_edges[ix] + eta_edges[ix + 1]) / 2.0),
                phi=float((phi_edges[iy] + phi_edges[iy + 1]) / 2.0),
                radius=float(np.hypot(eta_width / 2.0, phi_width / 2.0)),
                sigma=float(significance[ix, iy]),
            )
        )

    return FiducialMapSummary(
        mean_inefficiency=mean,
        stddev_inefficiency=stddev,
        inefficiency=inefficiency,
        significance=significance,
        hot_spots=tuple(hot_spots),
    )


def _walk_hists(value: Any):
    if hasattr(value, "axes") and hasattr(value, "values"):
        yield value
    elif isinstance(value, Mapping):
        for nested in value.values():
            yield from _walk_hists(nested)


def _hist2d_counts_edges(
    hist_obj: Any,
    *,
    category: str = "inclusive",
) -> tuple[np.ndarray, np.ndarray, np.ndarray] | None:
    try:
        axis_names = [axis.name for axis in hist_obj.axes]
    except AttributeError:
        return None

    if "cat" in axis_names:
        try:
            hist_obj = hist_obj[{"cat": category}]
        except Exception:
            return None

    try:
        axes = list(hist_obj.axes)
        counts = np.asarray(hist_obj.values(flow=False), dtype=float)
    except Exception:
        return None

    auxiliary_axis_names = {"cat", "variation", "sample"}
    keep_indices = [
        index
        for index, axis in enumerate(axes)
        if getattr(axis, "name", "") not in auxiliary_axis_names
        and hasattr(axis, "edges")
    ]
    if len(keep_indices) != 2:
        return None

    for index in reversed(range(len(axes))):
        if index not in keep_indices:
            counts = counts.sum(axis=index)

    eta_axis = axes[keep_indices[0]]
    phi_axis = axes[keep_indices[1]]
    eta_edges = np.asarray(eta_axis.edges, dtype=float)
    phi_edges = np.asarray(phi_axis.edges, dtype=float)
    expected_shape = (len(eta_edges) - 1, len(phi_edges) - 1)
    if counts.shape != expected_shape:
        return None
    return counts, eta_edges, phi_edges


def summed_hist2d_counts_edges(
    outputs: Sequence[Mapping[str, Any]],
    variable: str,
    *,
    dataset: str | None = None,
    sample: str | None = None,
    category: str = "inclusive",
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    total_counts = None
    total_eta_edges = None
    total_phi_edges = None
    for output in outputs:
        variables = output.get("variables", {})
        if variable not in variables:
            continue
        value: Any = variables[variable]
        if sample is not None and isinstance(value, Mapping):
            value = value.get(sample, {})
        if dataset is not None and sample is None and isinstance(value, Mapping):
            nested_values = value.values()
        elif dataset is not None and isinstance(value, Mapping):
            nested_values = [value.get(dataset, {})]
        else:
            nested_values = [value]

        for nested in nested_values:
            for hist_obj in _walk_hists(nested):
                result = _hist2d_counts_edges(hist_obj, category=category)
                if result is None:
                    continue
                counts, eta_edges, phi_edges = result
                if total_counts is None:
                    total_counts = counts.copy()
                    total_eta_edges = eta_edges.copy()
                    total_phi_edges = phi_edges.copy()
                    continue
                if (
                    len(eta_edges) != len(total_eta_edges)
                    or len(phi_edges) != len(total_phi_edges)
                    or not np.allclose(eta_edges, total_eta_edges)
                    or not np.allclose(phi_edges, total_phi_edges)
                ):
                    raise ValueError(f"histogram {variable!r} has inconsistent binning")
                total_counts += counts

    if total_counts is None or total_eta_edges is None or total_phi_edges is None:
        raise KeyError(f"2D histogram variable {variable!r} not found")
    return total_counts, total_eta_edges, total_phi_edges


def make_fiducial_map_from_outputs(
    outputs: Sequence[Mapping[str, Any]],
    *,
    before_variable: str,
    after_variable: str,
    dataset: str | None = None,
    sample: str | None = None,
    category: str = "inclusive",
    threshold: float = 2.0,
) -> tuple[FiducialMapSummary, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    before, eta_edges, phi_edges = summed_hist2d_counts_edges(
        outputs,
        before_variable,
        dataset=dataset,
        sample=sample,
        category=category,
    )
    after, after_eta_edges, after_phi_edges = summed_hist2d_counts_edges(
        outputs,
        after_variable,
        dataset=dataset,
        sample=sample,
        category=category,
    )
    if not np.allclose(eta_edges, after_eta_edges) or not np.allclose(
        phi_edges, after_phi_edges
    ):
        raise ValueError("before and after fiducial-map histograms have inconsistent binning")
    return (
        summarize_fiducial_map(before, after, eta_edges, phi_edges, threshold),
        before,
        after,
        eta_edges,
        phi_edges,
    )


def write_fiducial_map_payload(
    summary: FiducialMapSummary,
    *,
    before: np.ndarray,
    after: np.ndarray,
    eta_edges: np.ndarray,
    phi_edges: np.ndarray,
    output_json: Path,
    output_npz: Path | None = None,
    metadata: Mapping[str, Any] | None = None,
) -> Path:
    output_json.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "metadata": dict(metadata or {}),
        **summary.as_dict(),
    }
    output_json.write_text(json.dumps(payload, indent=2, sort_keys=True))

    if output_npz is not None:
        output_npz.parent.mkdir(parents=True, exist_ok=True)
        np.savez(
            output_npz,
            before=before,
            after=after,
            eta_edges=eta_edges,
            phi_edges=phi_edges,
            inefficiency=summary.inefficiency,
            significance=summary.significance,
        )
    return output_json


def _load_hot_spots(json_path: Path | None) -> tuple[HotSpot, ...]:
    if json_path is None:
        return ()
    payload = json.loads(json_path.read_text())
    return tuple(
        HotSpot(
            eta=float(hot_spot["eta"]),
            phi=float(hot_spot["phi"]),
            radius=float(hot_spot.get("radius", 0.06)),
            sigma=float(hot_spot["sigma"]),
        )
        for hot_spot in payload.get("hot_spots", [])
    )


def _z_range(flavor: str, quantity: str) -> tuple[float, float] | None:
    if quantity == "inefficiency":
        return (0.0, 0.5) if flavor == "electron" else (0.0, 0.05)
    if quantity == "significance_positive":
        return (0.0, 12.0) if flavor == "electron" else (0.0, 23.0)
    return None


def plot_fiducial_map_payload(
    npz_path: Path,
    *,
    output_prefix: Path,
    flavor: str,
    json_path: Path | None = None,
    run_period: str | None = None,
    lumi_text: str | None = None,
    cms_label: str = "CMS Preliminary",
    formats: Sequence[str] = ("pdf", "png"),
    draw_hot_spots: bool = True,
) -> list[Path]:
    try:
        import matplotlib.pyplot as plt
        from matplotlib.patches import Circle
    except ModuleNotFoundError as exc:
        raise RuntimeError("matplotlib is required to draw fiducial-map plots") from exc

    arrays = np.load(npz_path)
    eta_edges = arrays["eta_edges"]
    phi_edges = arrays["phi_edges"]
    hot_spots = _load_hot_spots(json_path) if draw_hot_spots else ()
    flavor_label = {"electron": "Electron", "muon": "Muon"}.get(flavor, flavor)
    period_suffix = "" if run_period is None else f" {run_period}"
    top_right = lumi_text or ""

    plots = (
        ("before", "beforeVeto", f"{flavor_label} fiducial map before veto{period_suffix}", "Events"),
        ("after", "afterVeto", f"{flavor_label} fiducial map after veto{period_suffix}", "Events"),
        (
            "inefficiency",
            "efficiency",
            f"{flavor_label} fiducial inefficiency{period_suffix}",
            "After / before",
        ),
        (
            "significance",
            "efficiencyInSigma",
            f"{flavor_label} fiducial inefficiency significance{period_suffix}",
            "Significance",
        ),
    )

    written = []
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    for key, legacy_name, title, colorbar_label in plots:
        values = np.asarray(arrays[key], dtype=float)
        if key == "significance":
            values = np.maximum(values, 0.0)
            zrange = _z_range(flavor, "significance_positive")
        else:
            zrange = _z_range(flavor, key)

        fig, ax = plt.subplots(figsize=(8, 8))
        mesh_kwargs = {"shading": "auto", "cmap": "viridis"}
        if zrange is not None:
            mesh_kwargs.update({"vmin": zrange[0], "vmax": zrange[1]})
        mesh = ax.pcolormesh(eta_edges, phi_edges, values.T, **mesh_kwargs)
        cbar = fig.colorbar(mesh, ax=ax, pad=0.015)
        cbar.set_label(colorbar_label, fontsize=13)
        cbar.ax.tick_params(labelsize=11)

        if key in ("inefficiency", "significance") and draw_hot_spots:
            for hot_spot in hot_spots:
                ax.add_patch(
                    Circle(
                        (hot_spot.eta, hot_spot.phi),
                        0.06,
                        fill=False,
                        linewidth=1.4,
                        edgecolor="#74c476",
                    )
                )

        ax.text(
            0.02,
            1.015,
            cms_label,
            transform=ax.transAxes,
            ha="left",
            va="bottom",
            fontsize=15,
            fontweight="bold",
        )
        if top_right:
            ax.text(
                0.98,
                1.015,
                top_right,
                transform=ax.transAxes,
                ha="right",
                va="bottom",
                fontsize=13,
            )
        ax.set_title(title, fontsize=14, pad=14)
        ax.set_xlabel(r"track $\eta$", fontsize=14)
        ax.set_ylabel(r"track $\phi$", fontsize=14)
        ax.set_xlim(float(eta_edges[0]), float(eta_edges[-1]))
        ax.set_ylim(float(phi_edges[0]), float(phi_edges[-1]))
        ax.tick_params(labelsize=12)
        fig.tight_layout()

        for fmt in formats:
            suffix = fmt if fmt.startswith(".") else f".{fmt}"
            out_path = output_prefix.parent / f"{output_prefix.name}_{legacy_name}{suffix}"
            fig.savefig(out_path, dpi=200)
            written.append(out_path)
        plt.close(fig)

    return written
