"""Small summary helpers for PocketCoffea output files."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from math import sqrt
from numbers import Number
from typing import Any


@dataclass(frozen=True)
class VetoProbabilitySummary:
    denominator_name: str
    numerator_name: str
    denominator: float
    numerator: float
    probability: float
    uncertainty: float

    def as_dict(self) -> dict[str, float | str]:
        return asdict(self)


def _sum_numeric_leaves(value: Any) -> float:
    """Sum all numeric leaves in a nested PocketCoffea cutflow object."""
    if isinstance(value, Number):
        return float(value)
    if isinstance(value, dict):
        return sum(_sum_numeric_leaves(v) for v in value.values())
    return 0.0


def cutflow_count(
    cutflow: dict[str, Any],
    category: str,
    *,
    dataset: str | None = None,
    sample: str | None = None,
    variation: str = "nominal",
) -> float:
    """Extract a category count from a nested PocketCoffea cutflow.

    PocketCoffea cutflows are nested differently depending on the stage.  The
    analysis categories we use here normally look like
    ``category -> dataset -> sample -> variation``.  This helper keeps the CLI
    forgiving: if filters are omitted it sums all matching numeric leaves.
    """
    if category not in cutflow:
        raise KeyError(f"category {category!r} not found in cutflow")

    value: Any = cutflow[category]
    if dataset is not None:
        value = value[dataset]
    if sample is not None:
        value = value[sample]
    if isinstance(value, dict) and variation in value:
        value = value[variation]
    return _sum_numeric_leaves(value)


def summarize_veto_probability(
    cutflow: dict[str, Any],
    *,
    denominator_name: str,
    numerator_name: str,
    dataset: str | None = None,
    sample: str | None = None,
    variation: str = "nominal",
) -> VetoProbabilitySummary:
    denominator = cutflow_count(
        cutflow,
        denominator_name,
        dataset=dataset,
        sample=sample,
        variation=variation,
    )
    numerator = cutflow_count(
        cutflow,
        numerator_name,
        dataset=dataset,
        sample=sample,
        variation=variation,
    )

    probability = numerator / denominator if denominator > 0.0 else 0.0
    uncertainty = (
        sqrt(probability * (1.0 - probability) / denominator)
        if denominator > 0.0
        else 0.0
    )

    return VetoProbabilitySummary(
        denominator_name=denominator_name,
        numerator_name=numerator_name,
        denominator=denominator,
        numerator=numerator,
        probability=probability,
        uncertainty=uncertainty,
    )
