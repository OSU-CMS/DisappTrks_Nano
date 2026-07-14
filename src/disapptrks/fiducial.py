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
