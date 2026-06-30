"""LaTeX table helpers for AN-style tag-and-probe summaries."""

from __future__ import annotations

from dataclasses import dataclass
from math import sqrt
from pathlib import Path
from typing import Any, Sequence

from .summaries import cutflow_count


@dataclass(frozen=True)
class CountWithVariance:
    value: float
    variance: float


@dataclass(frozen=True)
class AsymmetricVetoProbability:
    central: float
    err_down: float
    err_up: float
    numerator: float
    denominator: float


MUON_CUTFLOW_ROWS = [
    ("inclusive", r"Events processed"),
    ("muon_veto_tag", r"$\geq 1$ muon tag"),
    ("muon_veto_probe", r"$\geq 1$ muon tag and $\geq 1$ probe track"),
    ("muon_veto_pair", r"$\geq 1$ track--muon pair"),
    (
        "muon_veto_pair_os_mass10",
        r"$\geq 1$ opposite-sign track--muon pair with "
        r"$M_{\mathrm{track},\mu}>10~\mathrm{GeV}$",
    ),
    (
        "muon_veto_zwindow",
        r"$\geq 1$ opposite-sign track--muon pair with "
        r"$|M_{\mathrm{track},\mu}-M_Z|<10~\mathrm{GeV}$",
    ),
]

DISPLAY_LAYER = {
    "NLayers4": r"$N_{\mathrm{layers}}=4$",
    "NLayers5": r"$N_{\mathrm{layers}}=5$",
    "NLayers6plus": r"$N_{\mathrm{layers}}\geq 6$",
    "combinedBins": r"combined",
}

PVETO_TABLE_LAYERS = ("NLayers4", "NLayers5", "NLayers6plus", "combinedBins")


def format_count(value: float) -> str:
    if abs(value - round(value)) < 1.0e-9:
        return f"{round(value):d}"
    return f"{value:.3g}"


def format_pveto_latex(summary: AsymmetricVetoProbability) -> str:
    return (
        rf"${summary.central:.4g}^{{+{summary.err_up:.2g}}}"
        rf"_{{-{summary.err_down:.2g}}}$"
    )


def pveto_with_asymmetric_uncertainty(
    *,
    num_os: CountWithVariance,
    num_ss: CountWithVariance,
    den_os: CountWithVariance,
    den_ss: CountWithVariance,
) -> AsymmetricVetoProbability:
    """Compute the AN Pveto convention with SS subtraction.

    The convention follows the legacy table helper:

    ``P_veto = (N_veto_OS - N_veto_SS)/(N_OS - N_SS)``.

    If the signed numerator is negative, the central value is set to zero and
    only an upward one-sigma uncertainty is quoted.
    """
    numerator = num_os.value - num_ss.value
    denominator = den_os.value - den_ss.value

    if denominator <= 0.0:
        return AsymmetricVetoProbability(0.0, 0.0, 0.0, numerator, denominator)

    numerator_variance = max(num_os.variance + num_ss.variance, 0.0)
    denominator_variance = max(den_os.variance + den_ss.variance, 0.0)
    sigma_numerator = sqrt(numerator_variance)
    sigma_denominator = sqrt(denominator_variance)

    if numerator < 0.0:
        return AsymmetricVetoProbability(
            0.0,
            0.0,
            sigma_numerator / denominator,
            numerator,
            denominator,
        )

    central = numerator / denominator
    rel2 = 0.0
    if numerator > 0.0:
        rel2 += (sigma_numerator / numerator) ** 2
    if denominator > 0.0:
        rel2 += (sigma_denominator / denominator) ** 2
    uncertainty = abs(central) * sqrt(rel2)
    return AsymmetricVetoProbability(
        central,
        uncertainty,
        uncertainty,
        numerator,
        denominator,
    )


def _category_count(
    cutflow: dict[str, Any],
    category: str,
    *,
    dataset: str | None = None,
    sample: str | None = None,
    variation: str = "nominal",
) -> float:
    try:
        return cutflow_count(
            cutflow,
            category,
            dataset=dataset,
            sample=sample,
            variation=variation,
        )
    except KeyError:
        return 0.0


def write_muon_cutflow_latex(
    cutflow: dict[str, Any],
    path: Path,
    *,
    dataset: str | None = None,
    sample: str | None = None,
    variation: str = "nominal",
    include_table_env: bool = False,
) -> None:
    """Write an AN-style muon tag-and-probe cutflow table."""
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = [
        (
            label,
            _category_count(
                cutflow,
                category,
                dataset=dataset,
                sample=sample,
                variation=variation,
            ),
        )
        for category, label in MUON_CUTFLOW_ROWS
    ]

    with path.open("w") as out:
        if include_table_env:
            out.write(r"\begin{table}[htbp]" + "\n")
            out.write(r"\centering" + "\n")
            out.write(r"\caption{Muon tag-and-probe cutflow.}" + "\n")
            out.write(r"\label{tab:muon_tp_cutflow}" + "\n")

        out.write(r"\begin{tabular}{lrrr}" + "\n")
        out.write(r"\hline" + "\n")
        out.write(
            r"Cut & Events & $\epsilon_{\mathrm{prev}}$ & "
            r"$\epsilon_{\mathrm{total}}$ \\" + "\n"
        )
        out.write(r"\hline" + "\n")

        first = None
        previous = None
        for label, value in rows:
            if first is None:
                first = value
            eff_prev = value / previous if previous else 1.0
            eff_total = value / first if first else 0.0
            out.write(
                f"{label} & {format_count(value)} & "
                f"{eff_prev:.4f} & {eff_total:.4f} \\\\\n"
            )
            previous = value

        out.write(r"\hline" + "\n")
        out.write(r"\end{tabular}" + "\n")
        if include_table_env:
            out.write(r"\end{table}" + "\n")


def write_muon_pveto_latex(
    cutflow: dict[str, Any],
    path: Path,
    *,
    run_period: str,
    flavor: str = r"$\mu$",
    layers: Sequence[str] = PVETO_TABLE_LAYERS,
    dataset: str | None = None,
    sample: str | None = None,
    variation: str = "nominal",
    os_denominator_name: str = "muon_veto_zwindow",
    os_numerator_name: str = "muon_pveto_zwindow_pass",
    ss_denominator_name: str = "muon_veto_ss_zwindow",
    ss_numerator_name: str = "muon_pveto_ss_zwindow_pass",
    include_table_env: bool = False,
) -> dict[str, AsymmetricVetoProbability]:
    """Write AN/Table-24-style muon Pveto rows."""
    path.parent.mkdir(parents=True, exist_ok=True)

    if isinstance(layers, str):
        layers = (layers,)
    summaries = {}

    rows = []
    for layer in layers:
        suffix = "" if layer == "combinedBins" else f"_{layer}"
        den_os_value = _category_count(
            cutflow,
            f"{os_denominator_name}{suffix}",
            dataset=dataset,
            sample=sample,
            variation=variation,
        )
        num_os_value = _category_count(
            cutflow,
            f"{os_numerator_name}{suffix}",
            dataset=dataset,
            sample=sample,
            variation=variation,
        )
        den_ss_value = _category_count(
            cutflow,
            f"{ss_denominator_name}{suffix}",
            dataset=dataset,
            sample=sample,
            variation=variation,
        )
        num_ss_value = _category_count(
            cutflow,
            f"{ss_numerator_name}{suffix}",
            dataset=dataset,
            sample=sample,
            variation=variation,
        )

        summary = pveto_with_asymmetric_uncertainty(
            den_os=CountWithVariance(den_os_value, den_os_value),
            num_os=CountWithVariance(num_os_value, num_os_value),
            den_ss=CountWithVariance(den_ss_value, den_ss_value),
            num_ss=CountWithVariance(num_ss_value, num_ss_value),
        )
        summaries[layer] = summary
        rows.append((layer, den_os_value, num_os_value, den_ss_value, num_ss_value, summary))

    with path.open("w") as out:
        if include_table_env:
            out.write(r"\begin{table}[htbp]" + "\n")
            out.write(r"\centering" + "\n")
            out.write(r"\caption{Muon veto probability.}" + "\n")
            out.write(r"\label{tab:muon_pveto}" + "\n")

        out.write(r"\begin{tabular}{llcrrrrc}" + "\n")
        out.write(r"\hline" + "\n")
        out.write(
            r"run period & flavor & $n_{\mathrm{layers}}$ & "
            r"$N_{T\&P}$ & $N^{\mathrm{veto}}_{T\&P}$ & "
            r"$N_{SS,T\&P}$ & $N^{\mathrm{veto}}_{SS,T\&P}$ & "
            r"$P_{\mathrm{veto}}$ \\" + "\n"
        )
        out.write(r"\hline" + "\n")
        for layer, den_os_value, num_os_value, den_ss_value, num_ss_value, summary in rows:
            out.write(
                f"{run_period} & {flavor} & {DISPLAY_LAYER.get(layer, layer)} & "
                f"{format_count(den_os_value)} & {format_count(num_os_value)} & "
                f"{format_count(den_ss_value)} & {format_count(num_ss_value)} & "
                f"{format_pveto_latex(summary)} \\\\\n"
            )
        out.write(r"\hline" + "\n")
        out.write(r"\end{tabular}" + "\n")
        if include_table_env:
            out.write(r"\end{table}" + "\n")

    return summaries
