"""Fake-track background estimate helpers.

The Run-3 fake-track estimate follows the legacy DisappTrks convention

    N_fake = xi * N_ctrl

where ``xi`` is a d0 transfer factor measured from a 3-layer/basic-track d0
shape and ``N_ctrl`` is the yield in the target disappearing-track sideband.
When the control region is a Z->ll sample, ``N_ctrl`` is first normalized by the
ratio of the basic-search yield to the inclusive Z->ll yield.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from math import erf, sqrt
from pathlib import Path
from typing import Any, Mapping, Sequence

from .summaries import cutflow_count
from .tables import format_count, format_value_with_uncertainty


@dataclass(frozen=True)
class Count:
    value: float
    variance: float | None = None

    def __post_init__(self) -> None:
        if self.variance is None:
            object.__setattr__(self, "variance", max(float(self.value), 0.0))

    @property
    def error(self) -> float:
        return sqrt(max(float(self.variance), 0.0))

    def __mul__(self, other: "Count | float") -> "Count":
        if isinstance(other, Count):
            value = self.value * other.value
            variance = value * value * (_relative_variance(self) + _relative_variance(other))
            return Count(value, variance)
        return Count(self.value * other, float(self.variance) * other * other)

    def __truediv__(self, other: "Count | float") -> "Count":
        if isinstance(other, Count):
            if other.value == 0.0:
                return Count(0.0, 0.0)
            value = self.value / other.value
            variance = value * value * (_relative_variance(self) + _relative_variance(other))
            return Count(value, variance)
        if other == 0.0:
            return Count(0.0, 0.0)
        return Count(self.value / other, float(self.variance) / (other * other))


@dataclass(frozen=True)
class FakeTrackEstimate:
    layer: str
    transfer_signal: Count
    transfer_sideband: Count
    transfer_factor: Count
    control_raw: Count
    normalization: Count
    control: Count
    estimate: Count
    p_fake_raw: Count | None = None
    transfer_signal_category: str = ""
    transfer_sideband_category: str = ""
    control_category: str = ""
    basic_yield_category: str | None = None
    z_to_ll_yield_category: str | None = None

    def as_dict(self) -> dict[str, Any]:
        out = asdict(self)
        for key in (
            "transfer_signal",
            "transfer_sideband",
            "transfer_factor",
            "control_raw",
            "normalization",
            "control",
            "estimate",
            "p_fake_raw",
        ):
            value = out[key]
            if value is not None:
                out[key] = {
                    "value": value["value"],
                    "error": sqrt(max(value["variance"], 0.0)),
                    "variance": value["variance"],
                }
        return out


@dataclass(frozen=True)
class DxyTransferFactorFit:
    """Fit-derived d0 transfer factor used by the AN fake-track method."""

    control_region: str
    histogram: str
    numerator_range: tuple[float, float]
    denominator_range: tuple[float, float]
    fit_range: tuple[float, float]
    amplitude: float
    sigma: float
    constant: float
    transfer_factor: Count


@dataclass(frozen=True)
class ANFakeTrackEstimate:
    """Chapter-5 fake-track estimate components for one layer bin."""

    control_region: str
    layer: str
    control_events: Count
    sideband_events: Count
    basic_events: Count | None
    raw_probability: Count
    transfer_factor: Count
    fake_probability: Count
    fake_yield: Count | None
    control_category: str
    sideband_category: str
    basic_yield_category: str | None = None

    def as_dict(self) -> dict[str, Any]:
        out = asdict(self)
        for key in (
            "control_events",
            "sideband_events",
            "basic_events",
            "raw_probability",
            "transfer_factor",
            "fake_probability",
            "fake_yield",
        ):
            value = out[key]
            if value is not None:
                out[key] = {
                    "value": value["value"],
                    "error": sqrt(max(value["variance"], 0.0)),
                    "variance": value["variance"],
                }
        return out


def _relative_variance(count: Count) -> float:
    if count.value == 0.0:
        return 0.0
    return float(count.variance) / (count.value * count.value)


def _category_name(pattern: str | None, layer: str) -> str | None:
    if pattern is None:
        return None
    return pattern.format(layer=layer)


def _count_from_mapping(counts: Mapping[str, Any], category: str) -> Count:
    if category not in counts:
        raise KeyError(f"category {category!r} not found")
    value = counts[category]
    if isinstance(value, Mapping):
        if "value" not in value:
            raise KeyError(f"category {category!r} mapping must contain a 'value' key")
        return Count(float(value["value"]), float(value.get("variance", value["value"])))
    return Count(float(value))


def _count_from_cutflow(
    cutflow: dict[str, Any],
    category: str,
    *,
    dataset: str | None = None,
    sample: str | None = None,
    variation: str = "nominal",
) -> Count:
    value = cutflow_count(
        cutflow,
        category,
        dataset=dataset,
        sample=sample,
        variation=variation,
    )
    return Count(value)


def _find_hist_axis(hist_obj: Any):
    axes = list(getattr(hist_obj, "axes", []))
    for axis in axes:
        name = getattr(axis, "name", "")
        if name not in {"cat", "variation", "sample"}:
            return axis
    return axes[0] if axes else None


def _hist_counts_edges(hist_obj: Any, *, category: str = "inclusive") -> tuple[Any, Any] | None:
    """Return counts and bin edges from a one-dimensional PocketCoffea histogram."""

    try:
        axis_names = [axis.name for axis in hist_obj.axes]
    except AttributeError:
        return None

    if "cat" in axis_names:
        try:
            hist_obj = hist_obj[{"cat": category}]
        except Exception:
            return None

    axis = _find_hist_axis(hist_obj)
    if axis is None:
        return None

    try:
        import numpy as np

        values = hist_obj.values(flow=False)
        counts = np.asarray(values, dtype=float)
        edges = np.asarray(axis.edges, dtype=float)
    except Exception:
        return None

    while counts.ndim > 1:
        counts = counts.sum(axis=0)

    if counts.ndim != 1 or len(edges) != len(counts) + 1:
        return None
    return counts, edges


def _walk_hists(value: Any):
    if hasattr(value, "axes") and hasattr(value, "values"):
        yield value
    elif isinstance(value, Mapping):
        for nested in value.values():
            yield from _walk_hists(nested)


def summed_hist_counts_edges(
    outputs: Sequence[Mapping[str, Any]],
    variable: str,
    *,
    dataset: str | None = None,
    sample: str | None = None,
    category: str = "inclusive",
) -> tuple[Any, Any]:
    """Sum one histogram variable across PocketCoffea output files."""

    import numpy as np

    total_counts = None
    total_edges = None
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
                result = _hist_counts_edges(hist_obj, category=category)
                if result is None:
                    continue
                counts, edges = result
                if total_counts is None:
                    total_counts = counts.copy()
                    total_edges = edges.copy()
                else:
                    if len(edges) != len(total_edges) or not np.allclose(edges, total_edges):
                        raise ValueError(f"histogram {variable!r} has inconsistent binning")
                    total_counts += counts

    if total_counts is None or total_edges is None:
        raise KeyError(f"histogram variable {variable!r} not found")
    return total_counts, total_edges


def _gauss_plus_constant_integral(amplitude: float, sigma: float, constant: float, lo: float, hi: float) -> float:
    if sigma <= 0.0:
        return constant * (hi - lo)
    gaussian = amplitude * sigma * sqrt(3.141592653589793 / 2.0)
    gaussian *= erf(hi / (sqrt(2.0) * sigma)) - erf(lo / (sqrt(2.0) * sigma))
    return gaussian + constant * (hi - lo)


def _fit_gauss_plus_constant(x, y):
    import numpy as np
    from scipy.optimize import curve_fit

    def model(xvals, amplitude, sigma, constant):
        return amplitude * np.exp(-0.5 * (xvals / sigma) ** 2) + constant

    yerr = np.sqrt(np.maximum(y, 1.0))
    p0 = [max(float(y.max() - y.min()), 1.0), 0.10, max(float(y.min()), 0.0)]
    return curve_fit(
        model,
        x,
        y,
        p0=p0,
        sigma=yerr,
        absolute_sigma=True,
        bounds=([0.0, 1.0e-4, 0.0], [np.inf, 1.0, np.inf]),
        maxfev=20000,
    )


def _transfer_factor_from_fit(
    popt,
    pcov,
    *,
    numerator_ranges: Sequence[tuple[float, float]],
    denominator_ranges: Sequence[tuple[float, float]],
) -> Count:
    import numpy as np

    def ratio(params) -> float:
        numerator = sum(_gauss_plus_constant_integral(*params, *r) for r in numerator_ranges)
        denominator = sum(_gauss_plus_constant_integral(*params, *r) for r in denominator_ranges)
        return 0.0 if denominator == 0.0 else numerator / denominator

    value = ratio(popt)
    variance = 0.0
    try:
        grad = []
        for i, param in enumerate(popt):
            step = max(abs(float(param)) * 1.0e-5, 1.0e-7)
            shifted_hi = list(popt)
            shifted_lo = list(popt)
            shifted_hi[i] += step
            shifted_lo[i] -= step
            grad.append((ratio(shifted_hi) - ratio(shifted_lo)) / (2.0 * step))
        grad = np.asarray(grad)
        variance = float(grad @ pcov @ grad)
    except Exception:
        variance = 0.0
    return Count(value, max(variance, 0.0))


def fit_dxy_transfer_factor(
    counts,
    edges,
    *,
    control_region: str,
    histogram: str,
    numerator_range: tuple[float, float] = (0.0, 0.02),
    denominator_range: tuple[float, float] = (0.05, 0.50),
    fit_range: tuple[float, float] = (0.10, 0.50),
) -> DxyTransferFactorFit:
    """Fit |dxy| sideband with Gaussian(mean=0)+constant and return zeta."""

    import numpy as np

    centers = 0.5 * (edges[:-1] + edges[1:])
    fit_mask = (centers >= fit_range[0]) & (centers < fit_range[1])
    x = centers[fit_mask]
    y = counts[fit_mask]
    if len(x) < 3 or float(y.sum()) <= 0.0:
        raise ValueError(f"not enough entries to fit {histogram!r} for {control_region}")

    popt, pcov = _fit_gauss_plus_constant(x, y)
    amplitude, sigma, constant = (float(v) for v in popt)
    transfer_factor = _transfer_factor_from_fit(
        popt,
        pcov,
        numerator_ranges=[numerator_range],
        denominator_ranges=[denominator_range],
    )

    return DxyTransferFactorFit(
        control_region=control_region,
        histogram=histogram,
        numerator_range=numerator_range,
        denominator_range=denominator_range,
        fit_range=fit_range,
        amplitude=amplitude,
        sigma=sigma,
        constant=constant,
        transfer_factor=transfer_factor,
    )


def fit_signed_dxy_transfer_factor(
    counts,
    edges,
    *,
    control_region: str,
    histogram: str,
    numerator_abs_range: tuple[float, float] = (0.0, 0.02),
    denominator_abs_range: tuple[float, float] = (0.05, 0.50),
    fit_abs_range: tuple[float, float] = (0.10, 0.50),
) -> DxyTransferFactorFit:
    """Fit signed dxy sidebands, excluding central |dxy| < fit_abs_range[0]."""

    import numpy as np

    centers = 0.5 * (edges[:-1] + edges[1:])
    abs_centers = np.abs(centers)
    fit_mask = (abs_centers >= fit_abs_range[0]) & (abs_centers < fit_abs_range[1])
    x = centers[fit_mask]
    y = counts[fit_mask]
    if len(x) < 3 or float(y.sum()) <= 0.0:
        raise ValueError(f"not enough entries to fit {histogram!r} for {control_region}")

    popt, pcov = _fit_gauss_plus_constant(x, y)
    amplitude, sigma, constant = (float(v) for v in popt)
    numerator_ranges = [
        (-numerator_abs_range[1], -numerator_abs_range[0]),
        numerator_abs_range,
    ]
    denominator_ranges = [
        (-denominator_abs_range[1], -denominator_abs_range[0]),
        denominator_abs_range,
    ]
    transfer_factor = _transfer_factor_from_fit(
        popt,
        pcov,
        numerator_ranges=numerator_ranges,
        denominator_ranges=denominator_ranges,
    )

    return DxyTransferFactorFit(
        control_region=control_region,
        histogram=histogram,
        numerator_range=numerator_abs_range,
        denominator_range=denominator_abs_range,
        fit_range=fit_abs_range,
        amplitude=amplitude,
        sigma=sigma,
        constant=constant,
        transfer_factor=transfer_factor,
    )


def _dxy_fit_model(xvals, amplitude: float, sigma: float, constant: float):
    import numpy as np

    return amplitude * np.exp(-0.5 * (xvals / sigma) ** 2) + constant


def plot_dxy_transfer_factor(
    counts,
    edges,
    fit: DxyTransferFactorFit,
    path: Path,
    *,
    signed_counts=None,
    signed_edges=None,
    title: str | None = None,
) -> None:
    """Write a Figure-26-style signed-dxy histogram and transfer-factor fit plot.

    New outputs fit the signed ``dxy`` sidebands on both sides of zero while
    excluding the central fit gap.  Old outputs without signed histograms fall
    back to a mirrored visualization of the folded ``|dxy|`` histogram.
    """

    import numpy as np

    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise RuntimeError("matplotlib is required to write --fit-plot outputs") from exc

    path.parent.mkdir(parents=True, exist_ok=True)
    if signed_counts is None or signed_edges is None:
        # Backward-compatible fallback for old outputs.  The sign information is
        # not recoverable from |dxy|, so split the folded counts evenly into the
        # two signed halves.  New fake-track outputs should provide signed dxy.
        positive_edges = np.asarray(edges, dtype=float)
        negative_edges = -positive_edges[::-1]
        signed_edges = np.concatenate([negative_edges[:-1], positive_edges])
        signed_counts = 0.5 * np.concatenate([counts[::-1], counts])

    signed_counts = np.asarray(signed_counts, dtype=float)
    signed_edges = np.asarray(signed_edges, dtype=float)
    centers = 0.5 * (signed_edges[:-1] + signed_edges[1:])
    xfit = np.linspace(float(signed_edges[0]), float(signed_edges[-1]), 1000)
    yfit = _dxy_fit_model(np.abs(xfit), fit.amplitude, fit.sigma, fit.constant)

    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    ax.stairs(signed_counts, signed_edges, label="data", color="black", linewidth=1.5)
    ax.errorbar(
        centers,
        signed_counts,
        yerr=np.sqrt(np.maximum(signed_counts, 1.0)),
        fmt="o",
        color="black",
        markersize=3,
        linewidth=1,
    )
    ax.plot(
        xfit,
        yfit,
        color="tab:red",
        linewidth=2,
        label="Gaussian + constant fit",
    )

    for index, (lo, hi) in enumerate(
        ((-fit.numerator_range[1], -fit.numerator_range[0]), fit.numerator_range)
    ):
        ax.axvspan(
            lo,
            hi,
            color="tab:green",
            alpha=0.18,
            label="signal window" if index == 0 else None,
        )
    for index, (lo, hi) in enumerate(
        ((-fit.denominator_range[1], -fit.denominator_range[0]), fit.denominator_range)
    ):
        ax.axvspan(
            lo,
            hi,
            color="tab:blue",
            alpha=0.10,
            label="sideband" if index == 0 else None,
        )

    for boundary in (-fit.fit_range[1], -fit.fit_range[0], fit.fit_range[0], fit.fit_range[1]):
        ax.axvline(boundary, color="tab:red", linestyle="--", linewidth=1)

    ax.set_xlabel(r"track $d_{0}$ [cm]")
    ax.set_ylabel("Number of tracks / 0.04 cm")
    ax.set_xlim(float(signed_edges[0]), float(signed_edges[-1]))
    ax.set_ylim(bottom=0.0)
    ax.set_title(title or f"{fit.control_region} transfer-factor fit")
    ax.text(
        0.98,
        0.95,
        rf"$\zeta = {fit.transfer_factor.value:.3g} \pm {fit.transfer_factor.error:.2g}$"
        "\n"
        rf"$\sigma = {fit.sigma:.3g}$ cm",
        transform=ax.transAxes,
        ha="right",
        va="top",
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.8, "edgecolor": "0.8"},
    )
    ax.legend(loc="upper right", bbox_to_anchor=(1.0, 0.78), fontsize=8)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def estimate_fake_track_background_an(
    cutflow: Mapping[str, Any],
    *,
    layer: str,
    control_region: str,
    transfer_factor: Count,
    control_category: str,
    sideband_category: str,
    basic_yield_category: str | None = None,
    dataset: str | None = None,
    sample: str | None = None,
    variation: str = "nominal",
    basic_cutflow: Mapping[str, Any] | None = None,
    basic_dataset: str | None = None,
    basic_sample: str | None = None,
    basic_variation: str | None = None,
) -> ANFakeTrackEstimate:
    control_events = _count_from_cutflow(
        dict(cutflow),
        control_category,
        dataset=dataset,
        sample=sample,
        variation=variation,
    )
    sideband_events = _count_from_cutflow(
        dict(cutflow),
        sideband_category,
        dataset=dataset,
        sample=sample,
        variation=variation,
    )
    raw_probability = sideband_events / control_events
    fake_probability = raw_probability * transfer_factor

    basic_events = None
    fake_yield = None
    if basic_yield_category:
        basic_source = cutflow if basic_cutflow is None else basic_cutflow
        basic_events = _count_from_cutflow(
            dict(basic_source),
            basic_yield_category,
            dataset=dataset if basic_dataset is None else basic_dataset,
            sample=sample if basic_sample is None else basic_sample,
            variation=variation if basic_variation is None else basic_variation,
        )
        fake_yield = fake_probability * basic_events

    return ANFakeTrackEstimate(
        control_region=control_region,
        layer=layer,
        control_events=control_events,
        sideband_events=sideband_events,
        basic_events=basic_events,
        raw_probability=raw_probability,
        transfer_factor=transfer_factor,
        fake_probability=fake_probability,
        fake_yield=fake_yield,
        control_category=control_category,
        sideband_category=sideband_category,
        basic_yield_category=basic_yield_category,
    )


def estimate_fake_track_background(
    source: Mapping[str, Any],
    *,
    layer: str,
    transfer_signal_category: str = "fake_basic3hits_d0_signal",
    transfer_sideband_category: str = "fake_basic3hits_d0_sideband",
    control_category: str = "fake_control_{layer}",
    basic_yield_category: str | None = None,
    z_to_ll_yield_category: str | None = None,
    dataset: str | None = None,
    sample: str | None = None,
    variation: str = "nominal",
    source_is_cutflow: bool = True,
    prescale: float = 1.0,
) -> FakeTrackEstimate:
    """Compute the fake-track estimate for one layer bin.

    ``source`` can either be a PocketCoffea cutflow object or a flat mapping of
    category names to counts.  Category strings may include ``{layer}``, which is
    replaced by the current layer name.
    """

    def get_count(category_pattern: str) -> tuple[str, Count]:
        category = _category_name(category_pattern, layer)
        assert category is not None
        if source_is_cutflow:
            return category, _count_from_cutflow(
                source, category, dataset=dataset, sample=sample, variation=variation
            )
        return category, _count_from_mapping(source, category)

    transfer_signal_name, transfer_signal = get_count(transfer_signal_category)
    transfer_sideband_name, transfer_sideband = get_count(transfer_sideband_category)
    control_name, control_raw = get_count(control_category)

    transfer_factor = transfer_signal / transfer_sideband
    normalization = Count(1.0, 0.0)
    p_fake_raw = None
    basic_yield_name = _category_name(basic_yield_category, layer)
    z_to_ll_yield_name = _category_name(z_to_ll_yield_category, layer)

    if basic_yield_name is not None and z_to_ll_yield_name is not None:
        if source_is_cutflow:
            basic_yield = _count_from_cutflow(
                source,
                basic_yield_name,
                dataset=dataset,
                sample=sample,
                variation=variation,
            )
            z_to_ll_yield = _count_from_cutflow(
                source,
                z_to_ll_yield_name,
                dataset=dataset,
                sample=sample,
                variation=variation,
            )
        else:
            basic_yield = _count_from_mapping(source, basic_yield_name)
            z_to_ll_yield = _count_from_mapping(source, z_to_ll_yield_name)

        normalization = basic_yield / z_to_ll_yield
        p_fake_raw = control_raw / z_to_ll_yield

    control = control_raw * normalization * prescale
    estimate = control * transfer_factor

    return FakeTrackEstimate(
        layer=layer,
        transfer_signal=transfer_signal,
        transfer_sideband=transfer_sideband,
        transfer_factor=transfer_factor,
        control_raw=control_raw,
        normalization=normalization,
        control=control,
        estimate=estimate,
        p_fake_raw=p_fake_raw,
        transfer_signal_category=transfer_signal_name,
        transfer_sideband_category=transfer_sideband_name,
        control_category=control_name,
        basic_yield_category=basic_yield_name,
        z_to_ll_yield_category=z_to_ll_yield_name,
    )


def write_fake_track_latex(
    estimates: Sequence[FakeTrackEstimate],
    path: Path,
    *,
    run_period: str,
    include_table_env: bool = False,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as out:
        if include_table_env:
            out.write(r"\begin{table}[htbp]" + "\n")
            out.write(r"\centering" + "\n")
            out.write(r"\caption{Fake-track background estimate.}" + "\n")
            out.write(r"\label{tab:fake_track_estimate}" + "\n")

        out.write(r"\begin{tabular}{lrrrrr}" + "\n")
        out.write(r"\hline" + "\n")
        out.write(
            r"run period & $n_{\mathrm{layers}}$ & $N_{\mathrm{ctrl}}$ & "
            r"$\xi$ & $P_{\mathrm{fake}}^{\mathrm{raw}}$ & "
            r"$N_{\mathrm{fake}}$ \\" + "\n"
        )
        out.write(r"\hline" + "\n")
        for estimate in estimates:
            p_fake = (
                "--"
                if estimate.p_fake_raw is None
                else format_pm_latex(estimate.p_fake_raw.value, estimate.p_fake_raw.error)
            )
            out.write(
                f"{run_period} & {estimate.layer} & "
                f"{format_pm_latex(estimate.control.value, estimate.control.error)} & "
                f"{format_pm_latex(estimate.transfer_factor.value, estimate.transfer_factor.error)} & "
                f"{p_fake} & "
                f"{format_pm_latex(estimate.estimate.value, estimate.estimate.error)} \\\\\n"
            )
        out.write(r"\hline" + "\n")
        out.write(r"\end{tabular}" + "\n")
        if include_table_env:
            out.write(r"\end{table}" + "\n")


def write_an_fake_track_latex(
    estimates: Sequence[ANFakeTrackEstimate],
    fit: DxyTransferFactorFit,
    path: Path,
    *,
    run_period: str,
    include_table_env: bool = False,
) -> None:
    """Write the Chapter-5 fake-track estimate table."""

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as out:
        if include_table_env:
            out.write(r"\begin{table}[htbp]" + "\n")
            out.write(r"\centering" + "\n")
            out.write(r"\caption{Fake-track background estimate.}" + "\n")
            out.write(r"\label{tab:fake_track_estimate}" + "\n")

        has_yield = any(estimate.fake_yield is not None for estimate in estimates)
        if has_yield:
            out.write(r"\begin{tabular}{llrrrrrr}" + "\n")
            header = (
                r"run period & control & $n_{\mathrm{layers}}$ & $N_Z$ & "
                r"$N_{\mathrm{sideband}}$ & $\zeta$ & $P_{\mathrm{fake}}$ & "
                r"$N_{\mathrm{fake}}$ \\"
            )
        else:
            out.write(r"\begin{tabular}{llrrrrr}" + "\n")
            header = (
                r"run period & control & $n_{\mathrm{layers}}$ & $N_Z$ & "
                r"$N_{\mathrm{sideband}}$ & $\zeta$ & $P_{\mathrm{fake}}$ \\"
            )
        out.write(r"\hline" + "\n")
        out.write(header + "\n")
        out.write(r"\hline" + "\n")

        for estimate in estimates:
            zeta = format_pm_latex(
                estimate.transfer_factor.value,
                estimate.transfer_factor.error,
            )
            p_fake = format_pm_latex(
                estimate.fake_probability.value,
                estimate.fake_probability.error,
            )
            row = (
                f"{run_period} & {estimate.control_region} & {estimate.layer} & "
                f"{format_count(estimate.control_events.value)} & "
                f"{format_count(estimate.sideband_events.value)} & "
                f"{zeta} & {p_fake}"
            )
            if has_yield:
                if estimate.fake_yield is None:
                    row += " & --"
                else:
                    row += f" & {format_pm_latex(estimate.fake_yield.value, estimate.fake_yield.error)}"
            out.write(row + r" \\" + "\n")

        out.write(r"\hline" + "\n")
        out.write(r"\end{tabular}" + "\n")
        out.write(
            "% Fit: Gaussian with mean fixed to zero plus constant, "
            f"{fit.fit_range[0]} <= |dxy| < {fit.fit_range[1]} cm. "
            f"sigma={fit.sigma:.4g}, constant={fit.constant:.4g}.\n"
        )
        if include_table_env:
            out.write(r"\end{table}" + "\n")


def _count_from_payload(value: Mapping[str, Any]) -> Count:
    return Count(float(value["value"]), float(value.get("variance", value.get("error", 0.0) ** 2)))


def _layer_label(layer: str) -> str:
    return {
        "NLayers4": "4",
        "NLayers5": "5",
        "NLayers6plus": r"$\geq 6$",
        "combinedBins": "combined",
    }.get(layer, layer)


def _control_column_key(control_region: str) -> str:
    lower = control_region.lower()
    if r"\mu" in lower or "mumu" in lower or "mu" in lower:
        return "zmumu"
    if "ee" in lower or "electron" in lower:
        return "zee"
    return control_region


def _format_scientific_pm(value: float, error: float) -> str:
    if value == 0.0 and error == 0.0:
        return r"$0 \pm 0$"

    import math

    scale_value = abs(value) if value != 0.0 else abs(error)
    exponent = int(math.floor(math.log10(scale_value))) if scale_value > 0.0 else 0
    mantissa = value / (10.0**exponent)
    mantissa_error = error / (10.0**exponent)
    mantissa_text, error_text = format_value_with_uncertainty(
        mantissa,
        mantissa_error,
    )
    return rf"$({mantissa_text} \pm {error_text}) \times 10^{{{exponent}}}$"


def _format_yield_pm(value: float, error: float) -> str:
    value_text, error_text = format_value_with_uncertainty(value, error)
    return rf"${value_text} \pm {error_text}$"


def write_fake_track_table34_latex(
    json_paths: Sequence[Path],
    path: Path,
    *,
    run_period: str,
    include_table_env: bool = False,
) -> None:
    """Write an AN Table-34-style comparison of Z->mumu and Z->ee estimates."""

    by_layer: dict[str, dict[str, dict[str, Count]]] = {}
    for json_path in json_paths:
        payload = json.loads(json_path.read_text())
        for estimate in payload.get("estimates", []):
            layer = estimate["layer"]
            control = _control_column_key(estimate["control_region"])
            by_layer.setdefault(layer, {})[control] = {
                "p_fake": _count_from_payload(estimate["fake_probability"]),
                "n_fake": _count_from_payload(estimate["fake_yield"])
                if estimate.get("fake_yield") is not None
                else Count(0.0, 0.0),
            }

    layer_order = ["NLayers4", "NLayers5", "NLayers6plus", "combinedBins"]
    layers = [layer for layer in layer_order if layer in by_layer] + [
        layer for layer in by_layer if layer not in layer_order
    ]

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as out:
        if include_table_env:
            out.write(r"\begin{table}[htbp]" + "\n")
            out.write(r"\centering" + "\n")
            out.write(r"\caption{Comparison of fake-track event probabilities and background estimates.}" + "\n")
            out.write(r"\label{tab:fake_track_table34}" + "\n")

        out.write(r"\begin{tabular}{llrrrr}" + "\n")
        out.write(r"\hline" + "\n")
        out.write(
            r"run period & $n_{\mathrm{layers}}$ & "
            r"$P_{\mathrm{fake}}(Z\to\mu\mu)$ & "
            r"$P_{\mathrm{fake}}(Z\to ee)$ & "
            r"$N^{\mathrm{fake}}_{\mathrm{est}}(Z\to\mu\mu)$ & "
            r"$N^{\mathrm{fake}}_{\mathrm{est}}(Z\to ee)$ \\" + "\n"
        )
        out.write(r"\hline" + "\n")
        for i, layer in enumerate(layers):
            controls = by_layer[layer]
            zmumu = controls.get("zmumu")
            zee = controls.get("zee")

            def get(control: dict[str, Count] | None, key: str, formatter) -> str:
                if control is None:
                    return "--"
                count = control[key]
                return formatter(count.value, count.error)

            out.write(
                f"{run_period if i == 0 else ''} & {_layer_label(layer)} & "
                f"{get(zmumu, 'p_fake', _format_scientific_pm)} & "
                f"{get(zee, 'p_fake', _format_scientific_pm)} & "
                f"{get(zmumu, 'n_fake', _format_yield_pm)} & "
                f"{get(zee, 'n_fake', _format_yield_pm)} \\\\\n"
            )
        out.write(r"\hline" + "\n")
        out.write(r"\end{tabular}" + "\n")
        if include_table_env:
            out.write(r"\end{table}" + "\n")
