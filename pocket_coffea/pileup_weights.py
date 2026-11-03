"""Pileup weight reweighting Summer23BPix (2023D) MC to a later year's data.

The 2024-2026 signal yields use 2023D signal MC, and the missing-hits
correction fits for those years compare 2023D background MC with that year's
data.  PocketCoffea's built-in ``pileup`` weight looks up the LUM payload by
the sample's own year, i.e. reweights to 2023D data.  ``pileup_target``
instead reweights Summer23BPix MC to the target year's data pileup profile,
using the payload built by ``scripts/make_target_pileup_weights.py`` (2026 is
a labelled stand-in using the 2025 profile until LUM publishes one).

The payload is read once at configuration time and embedded in the weight,
so Dask workers need no extra transferred file.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import correctionlib

from pocket_coffea.lib.weights.weights import WeightLambda

PAYLOAD_PATH = (
    Path(__file__).resolve().parent
    / "data"
    / "pileup"
    / "puWeights_Summer23BPix_to_target.json.gz"
)
SOURCE_YEAR = "2023_postBPix"
TARGET_YEARS = ("2024", "2025", "2026")


def load_payload(path: Path = PAYLOAD_PATH) -> str:
    with gzip.open(path, "rt") as handle:
        return handle.read()


def make_pileup_target_weight(target_year: str, payload: str):
    """Return the ``pileup_target`` weight class for ``target_year``."""

    if target_year not in TARGET_YEARS:
        raise ValueError(
            f"DISAPPTRKS_PILEUP_TARGET_YEAR={target_year!r}: expected one of "
            f"{', '.join(TARGET_YEARS)}"
        )
    correction_name = f"Summer23BPix_to_{target_year}"

    def compute(params, metadata, events, size, shape_variations):
        if metadata["year"] != SOURCE_YEAR:
            raise ValueError(
                f"pileup_target ({correction_name}) is only valid for "
                f"{SOURCE_YEAR} (Summer23BPix) MC, got year={metadata['year']!r}"
            )
        correction = correctionlib.CorrectionSet.from_string(payload)[correction_name]
        n_true = events.Pileup.nTrueInt.to_numpy()
        return (
            correction.evaluate(n_true, "nominal"),
            correction.evaluate(n_true, "up"),
            correction.evaluate(n_true, "down"),
        )

    return WeightLambda.wrap_func(
        name="pileup_target",
        function=compute,
        has_variations=True,
    )
