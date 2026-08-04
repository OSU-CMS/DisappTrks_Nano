"""Lepton-background estimates built from Pveto-style control counts."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from math import sqrt
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

from .fake_tracks import Count
from .summaries import cutflow_count
from .tables import (
    CountWithVariance,
    POISSON_ZERO_UPPER_68,
    format_pm_latex,
    pveto_with_asymmetric_uncertainty,
)


@dataclass(frozen=True)
class LeptonBackgroundEstimate:
    """One layer-bin lepton-background estimate."""

    flavor: str
    layer: str
    control_raw: Count
    control: Count
    p_veto: Count
    p_offline: Count
    p_miss: Count
    trigger_efficiency: Count
    estimate: Count
    control_category: str
    poffline_numerator_category: str
    poffline_denominator_category: str
    pmiss_numerator_category: str
    pmiss_denominator_category: str

    def as_dict(self) -> dict[str, Any]:
        out = asdict(self)
        for key in (
            "control_raw",
            "control",
            "p_veto",
            "p_offline",
            "p_miss",
            "trigger_efficiency",
            "estimate",
        ):
            value = out[key]
            out[key] = {
                "value": value["value"],
                "error": sqrt(max(value["variance"], 0.0)),
                "variance": value["variance"],
            }
        return out

    @property
    def p_trigger(self) -> Count:
        """Backward-compatible alias for the AN ``Pmiss`` factor."""

        return self.p_miss


def _category_name(pattern: str, layer: str) -> str:
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


def _count_from_source(
    *,
    counts: Mapping[str, Any] | None,
    cutflow: dict[str, Any] | None,
    category: str,
    dataset: str | None,
    sample: str | None,
    variation: str,
) -> Count:
    if counts is not None:
        return _count_from_mapping(counts, category)
    if cutflow is None:
        raise ValueError("either counts or cutflow must be provided")
    return _count_from_cutflow(
        cutflow,
        category,
        dataset=dataset,
        sample=sample,
        variation=variation,
    )


def probability_from_counts(numerator: Count, denominator: Count) -> Count:
    """Return a binomial probability with a conservative zero-denominator guard."""

    if denominator.value <= 0.0:
        return Count(0.0, 0.0)
    value = numerator.value / denominator.value
    clipped = min(max(value, 0.0), 1.0)
    variance = clipped * (1.0 - clipped) / denominator.value
    return Count(value, variance)


def trigger_efficiency_from_counts(
    *,
    total_os: Count,
    total_ss: Count,
    passes_os: Count,
    passes_ss: Count,
) -> Count:
    """Legacy OS-minus-SS trigger-efficiency ratio helper."""

    passes = Count(
        passes_os.value - passes_ss.value,
        passes_os.variance + passes_ss.variance,
    )
    total = Count(
        total_os.value - total_ss.value,
        total_os.variance + total_ss.variance,
    )
    return probability_from_counts(passes, total)


def _walk_hists(value: Any):
    if hasattr(value, "axes") and hasattr(value, "values"):
        yield value
    elif isinstance(value, Mapping):
        for nested in value.values():
            yield from _walk_hists(nested)


def _iter_variable_hists(
    outputs: Sequence[Mapping[str, Any]],
    variable: str,
    *,
    dataset: str | None = None,
    sample: str | None = None,
):
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
            yield from _walk_hists(nested)


def _strip_category_axes(hist_obj: Any, *, category: str = "inclusive") -> Any | None:
    try:
        axis_names = [axis.name for axis in hist_obj.axes]
    except AttributeError:
        return None

    for axis_name in ("cat", "variation", "sample"):
        if axis_name in axis_names:
            try:
                selector = category if axis_name == "cat" else "nominal"
                hist_obj = hist_obj[{axis_name: selector}]
            except Exception:
                if axis_name == "cat":
                    return None
    return hist_obj


def _hist_counts_edges_nd(hist_obj: Any, *, category: str = "inclusive"):
    hist_obj = _strip_category_axes(hist_obj, category=category)
    if hist_obj is None:
        return None
    axes = [axis for axis in hist_obj.axes if getattr(axis, "name", "") != "cat"]
    try:
        counts = np.asarray(hist_obj.values(flow=False), dtype=float)
        edges = [np.asarray(axis.edges, dtype=float) for axis in axes]
    except Exception:
        return None
    if counts.ndim != len(edges):
        return None
    return counts, edges


def _summed_hist_counts_edges_nd(
    outputs: Sequence[Mapping[str, Any]],
    variable: str,
    *,
    dataset: str | None = None,
    sample: str | None = None,
    category: str = "inclusive",
):
    total_counts = None
    total_edges = None
    for hist_obj in _iter_variable_hists(
        outputs,
        variable,
        dataset=dataset,
        sample=sample,
    ):
        result = _hist_counts_edges_nd(hist_obj, category=category)
        if result is None:
            continue
        counts, edges = result
        if total_counts is None:
            total_counts = counts.copy()
            total_edges = [edge.copy() for edge in edges]
        else:
            if counts.shape != total_counts.shape:
                raise ValueError(f"incompatible histogram shapes for {variable}")
            total_counts += counts
    if total_counts is None or total_edges is None:
        return None
    return total_counts, total_edges


def _hist_integral(
    counts: np.ndarray,
    edges: Sequence[np.ndarray],
    *,
    met_cut: float,
    phi_cut: float,
) -> float:
    if counts.ndim != 2 or len(edges) != 2:
        return 0.0
    met_centers = 0.5 * (edges[0][:-1] + edges[0][1:])
    phi_centers = 0.5 * (edges[1][:-1] + edges[1][1:])
    mask = (met_centers[:, None] >= met_cut) & (phi_centers[None, :] >= phi_cut)
    return float(counts[mask].sum())


def _weighted_met_trigger_integral(
    met_phi_counts: np.ndarray,
    met_phi_edges: Sequence[np.ndarray],
    trigger_total_counts: np.ndarray,
    trigger_pass_counts: np.ndarray,
    trigger_edges: Sequence[np.ndarray],
    *,
    met_cut: float,
    phi_cut: float,
) -> tuple[float, float]:
    if met_phi_counts.ndim != 2 or len(met_phi_edges) != 2:
        return 0.0, 0.0
    if trigger_total_counts.ndim != 1 or trigger_pass_counts.ndim != 1:
        return 0.0, 0.0

    met_centers = 0.5 * (met_phi_edges[0][:-1] + met_phi_edges[0][1:])
    phi_centers = 0.5 * (met_phi_edges[1][:-1] + met_phi_edges[1][1:])
    trigger_centers = 0.5 * (trigger_edges[0][:-1] + trigger_edges[0][1:])
    with np.errstate(divide="ignore", invalid="ignore"):
        turn_on = np.divide(
            trigger_pass_counts,
            trigger_total_counts,
            out=np.zeros_like(trigger_pass_counts, dtype=float),
            where=trigger_total_counts > 0.0,
        )
    met_bin_turn_on = np.interp(
        met_centers,
        trigger_centers,
        turn_on,
        left=0.0,
        right=turn_on[-1] if len(turn_on) else 0.0,
    )
    offline_mask = (met_centers[:, None] >= met_cut) & (
        phi_centers[None, :] >= phi_cut
    )
    offline_total = float(met_phi_counts[offline_mask].sum())
    weighted = float((met_phi_counts * met_bin_turn_on[:, None])[offline_mask].sum())
    return weighted, offline_total


def legacy_met_probabilities_from_outputs(
    outputs: Sequence[Mapping[str, Any]],
    *,
    prefix: str,
    layers: Sequence[str],
    control_counts: Mapping[str, Count],
    dataset: str | None = None,
    sample: str | None = None,
    category: str = "inclusive",
    met_cut: float = 120.0,
    phi_cut: float = 0.5,
) -> dict[str, tuple[Count, Count]]:
    """Compute Poffline and Pmiss with the legacy MET-trigger turn-on method."""

    probabilities: dict[str, tuple[Count, Count]] = {}
    for layer in layers:
        met_variable = f"n{prefix}BackgroundMetNoMuPt_{layer}"
        trig_variable = f"n{prefix}BackgroundMetNoMuPtTrig_{layer}"
        fallback_met_variable = f"n{prefix}BackgroundMetMinusOnePt_{layer}"
        fallback_trig_variable = f"n{prefix}BackgroundMetMinusOnePtTrig_{layer}"
        met_phi_variable = (
            f"n{prefix}BackgroundDeltaPhiMetJetLeadingVsMetMinusOnePt_{layer}"
        )
        met_hist = _summed_hist_counts_edges_nd(
            outputs,
            met_variable,
            dataset=dataset,
            sample=sample,
            category=category,
        )
        trig_hist = _summed_hist_counts_edges_nd(
            outputs,
            trig_variable,
            dataset=dataset,
            sample=sample,
            category=category,
        )
        if met_hist is None or trig_hist is None:
            met_hist = _summed_hist_counts_edges_nd(
                outputs,
                fallback_met_variable,
                dataset=dataset,
                sample=sample,
                category=category,
            )
            trig_hist = _summed_hist_counts_edges_nd(
                outputs,
                fallback_trig_variable,
                dataset=dataset,
                sample=sample,
                category=category,
            )
        met_phi_hist = _summed_hist_counts_edges_nd(
            outputs,
            met_phi_variable,
            dataset=dataset,
            sample=sample,
            category=category,
        )
        if met_hist is None or trig_hist is None or met_phi_hist is None:
            continue
        met_counts, met_edges = met_hist
        trig_counts, trig_edges = trig_hist
        met_phi_counts, met_phi_edges = met_phi_hist

        offline_pass = _hist_integral(
            met_phi_counts,
            met_phi_edges,
            met_cut=met_cut,
            phi_cut=phi_cut,
        )
        poffline = probability_from_counts(
            Count(offline_pass, offline_pass),
            control_counts[layer],
        )
        weighted_trigger_pass, offline_total = _weighted_met_trigger_integral(
            met_phi_counts,
            met_phi_edges,
            met_counts,
            trig_counts,
            met_edges,
            met_cut=met_cut,
            phi_cut=phi_cut,
        )
        pmiss = probability_from_counts(
            Count(weighted_trigger_pass, weighted_trigger_pass),
            Count(offline_total, offline_total),
        )
        probabilities[layer] = (poffline, pmiss)
    return probabilities


def pveto_count_from_pair_counts(
    pair_counts: Mapping[str, float],
    *,
    use_two_lepton_denominator: bool = False,
) -> Count:
    """Convert OS/SS Pveto pair counts into a symmetric Count approximation.

    The legacy histogram-based path, used by the 2022/2023 background scripts,
    uses the same-sign-subtracted direct ratio.  The two-lepton denominator is
    kept for reproducing the older non-histogram fallback path.
    """

    if use_two_lepton_denominator:
        numerator = pair_counts.get("num_os", 0.0) - pair_counts.get("num_ss", 0.0)
        denominator = (
            2.0 * (pair_counts.get("den_os", 0.0) - pair_counts.get("den_ss", 0.0))
            - numerator
        )
        if denominator <= 0.0:
            return Count(0.0, 0.0)
        if numerator <= 0.0:
            return Count(0.0, (POISSON_ZERO_UPPER_68 / denominator) ** 2)
        value = numerator / denominator
        variance = value * (1.0 - min(max(value, 0.0), 1.0)) / denominator
        return Count(value, variance)

    summary = pveto_with_asymmetric_uncertainty(
        den_os=CountWithVariance(pair_counts.get("den_os", 0.0), pair_counts.get("den_os", 0.0)),
        num_os=CountWithVariance(pair_counts.get("num_os", 0.0), pair_counts.get("num_os", 0.0)),
        den_ss=CountWithVariance(pair_counts.get("den_ss", 0.0), pair_counts.get("den_ss", 0.0)),
        num_ss=CountWithVariance(pair_counts.get("num_ss", 0.0), pair_counts.get("num_ss", 0.0)),
    )
    uncertainty = max(summary.err_down, summary.err_up)
    return Count(summary.central, uncertainty * uncertainty)


def estimate_lepton_background(
    *,
    flavor: str,
    layers: Sequence[str],
    pair_counts: Mapping[str, Mapping[str, float]],
    cutflow: dict[str, Any] | None = None,
    counts: Mapping[str, Any] | None = None,
    control_category: str,
    poffline_numerator_category: str,
    poffline_denominator_category: str,
    pmiss_numerator_category: str | None = None,
    pmiss_denominator_category: str | None = None,
    ptrigger_numerator_category: str | None = None,
    ptrigger_denominator_category: str | None = None,
    control_prescale: float = 1.0,
    trigger_efficiency: Count | Mapping[str, Count] | None = None,
    met_probabilities: Mapping[str, tuple[Count, Count]] | None = None,
    dataset: str | None = None,
    sample: str | None = None,
    variation: str = "nominal",
) -> list[LeptonBackgroundEstimate]:
    estimates = []
    pmiss_numerator_category = pmiss_numerator_category or ptrigger_numerator_category
    pmiss_denominator_category = pmiss_denominator_category or ptrigger_denominator_category
    if pmiss_numerator_category is None or pmiss_denominator_category is None:
        raise ValueError(
            "pmiss_numerator_category and pmiss_denominator_category are required"
        )
    trigger_efficiency = trigger_efficiency or Count(1.0, 0.0)
    for layer in layers:
        control_name = _category_name(control_category, layer)
        poffline_num_name = _category_name(poffline_numerator_category, layer)
        poffline_den_name = _category_name(poffline_denominator_category, layer)
        pmiss_num_name = _category_name(pmiss_numerator_category, layer)
        pmiss_den_name = _category_name(pmiss_denominator_category, layer)

        control_raw = _count_from_source(
            counts=counts,
            cutflow=cutflow,
            category=control_name,
            dataset=dataset,
            sample=sample,
            variation=variation,
        )
        control = control_raw * control_prescale
        if met_probabilities is not None and layer in met_probabilities:
            poffline, pmiss = met_probabilities[layer]
        else:
            poffline = probability_from_counts(
                _count_from_source(
                    counts=counts,
                    cutflow=cutflow,
                    category=poffline_num_name,
                    dataset=dataset,
                    sample=sample,
                    variation=variation,
                ),
                _count_from_source(
                    counts=counts,
                    cutflow=cutflow,
                    category=poffline_den_name,
                    dataset=dataset,
                    sample=sample,
                    variation=variation,
                ),
            )
            pmiss = probability_from_counts(
                _count_from_source(
                    counts=counts,
                    cutflow=cutflow,
                    category=pmiss_num_name,
                    dataset=dataset,
                    sample=sample,
                    variation=variation,
                ),
                _count_from_source(
                    counts=counts,
                    cutflow=cutflow,
                    category=pmiss_den_name,
                    dataset=dataset,
                    sample=sample,
                    variation=variation,
                ),
            )
        p_veto = pveto_count_from_pair_counts(
            pair_counts.get(layer, {}),
            use_two_lepton_denominator=False,
        )
        layer_trigger_efficiency = (
            trigger_efficiency.get(layer, Count(1.0, 0.0))
            if isinstance(trigger_efficiency, Mapping)
            else trigger_efficiency
        )
        estimates.append(
            LeptonBackgroundEstimate(
                flavor=flavor,
                layer=layer,
                control_raw=control_raw,
                control=control,
                p_veto=p_veto,
                p_offline=poffline,
                p_miss=pmiss,
                trigger_efficiency=layer_trigger_efficiency,
                estimate=control
                * p_veto
                * poffline
                * pmiss
                / layer_trigger_efficiency,
                control_category=control_name,
                poffline_numerator_category=poffline_num_name,
                poffline_denominator_category=poffline_den_name,
                pmiss_numerator_category=pmiss_num_name,
                pmiss_denominator_category=pmiss_den_name,
            )
        )
    return estimates


def write_lepton_background_json(
    estimates: Sequence[LeptonBackgroundEstimate],
    path: Path,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps([estimate.as_dict() for estimate in estimates], indent=2, sort_keys=True))


def _an_layer_label(layer: str) -> str:
    return {
        "NLayers4": "4",
        "NLayers5": "5",
        "NLayers6plus": r"$\geq 6$",
        "combinedBins": "combined",
    }.get(layer, layer)


def _an_background_title(flavor: str) -> str:
    normalized = flavor.replace("$", "").replace("\\", "").lower()
    if "tau" in normalized:
        return "tau background"
    if "mu" in normalized:
        return "muon background"
    if "e" in normalized:
        return "electron background"
    return f"{flavor} background"


def write_lepton_background_latex(
    estimates: Sequence[LeptonBackgroundEstimate],
    path: Path,
    *,
    run_period: str,
    include_table_env: bool = False,
    tau_probability: Count | None = None,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    title = _an_background_title(estimates[0].flavor) if estimates else "lepton background"
    include_tau_probability = tau_probability is not None
    n_columns = 9 if include_tau_probability else 8
    column_spec = "ll" + ("c" * (n_columns - 2))
    with path.open("w") as out:
        if include_table_env:
            out.write(r"\begin{table}" + "\n")
            out.write(r"\centering" + "\n")
            out.write(r"\caption{" + title.capitalize() + r" estimate.}" + "\n")
            out.write(r"\label{tab:lepton_background_estimate}" + "\n")
        out.write(r"\begin{tiny}" + "\n")
        out.write(r"\begin{tabular}{" + column_spec + r"}" + "\n")
        out.write(r"\hline" + "\n")
        out.write(r"\multicolumn{" + str(n_columns) + r"}{c}{" + title + r"} \\" + "\n")
        out.write(
            r"run period & $n_{\mathrm{layers}}$ & "
            r"$\epsilon_{\mathrm{trigger}}^{\ell}$ & "
            + (r"$P(\tau)$ & " if include_tau_probability else "")
            + r"$N_{\mathrm{ctrl}}^{\ell}$ & $P_{\mathrm{veto}}$ & "
            r"$P_{\mathrm{offline}}$ & $P_{\mathrm{trigger}}$ & estimate \\" + "\n"
        )
        out.write(r"\hline" + "\n")
        n_rows = len(estimates)
        for index, estimate in enumerate(estimates):
            run_period_cell = rf"\multirow{{{n_rows}}}{{*}}{{{run_period}}}" if index == 0 else ""
            tau_probability_cell = (
                f"{format_pm_latex(tau_probability.value, tau_probability.error)} & "
                if include_tau_probability
                else ""
            )
            out.write(
                f"{run_period_cell} & {_an_layer_label(estimate.layer)} & "
                f"{format_pm_latex(estimate.trigger_efficiency.value, estimate.trigger_efficiency.error)} & "
                f"{tau_probability_cell}"
                f"{format_pm_latex(estimate.control.value, estimate.control.error)} & "
                f"{format_pm_latex(estimate.p_veto.value, estimate.p_veto.error)} & "
                f"{format_pm_latex(estimate.p_offline.value, estimate.p_offline.error)} & "
                f"{format_pm_latex(estimate.p_miss.value, estimate.p_miss.error)} & "
                f"{format_pm_latex(estimate.estimate.value, estimate.estimate.error)} "
                r"\\" + "\n"
            )
        out.write(r"\hline" + "\n")
        out.write(r"\end{tabular}" + "\n")
        out.write(r"\end{tiny}" + "\n")
        if include_table_env:
            out.write(r"\end{table}" + "\n")
