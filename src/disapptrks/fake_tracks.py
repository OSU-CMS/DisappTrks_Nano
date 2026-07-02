"""Fake-track background estimate helpers.

The Run-3 fake-track estimate follows the legacy DisappTrks convention

    N_fake = xi * N_ctrl

where ``xi`` is a d0 transfer factor measured from a 3-layer/basic-track d0
shape and ``N_ctrl`` is the yield in the target disappearing-track sideband.
When the control region is a Z->ll sample, ``N_ctrl`` is first normalized by the
ratio of the basic-search yield to the inclusive Z->ll yield.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from math import sqrt
from pathlib import Path
from typing import Any, Mapping, Sequence

from .summaries import cutflow_count
from .tables import format_count


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
                else f"{estimate.p_fake_raw.value:.4g} $\\pm$ {estimate.p_fake_raw.error:.2g}"
            )
            out.write(
                f"{run_period} & {estimate.layer} & "
                f"{format_count(estimate.control.value)} $\\pm$ {estimate.control.error:.2g} & "
                f"{estimate.transfer_factor.value:.4g} $\\pm$ {estimate.transfer_factor.error:.2g} & "
                f"{p_fake} & "
                f"{format_count(estimate.estimate.value)} $\\pm$ {estimate.estimate.error:.2g} \\\\\n"
            )
        out.write(r"\hline" + "\n")
        out.write(r"\end{tabular}" + "\n")
        if include_table_env:
            out.write(r"\end{table}" + "\n")
