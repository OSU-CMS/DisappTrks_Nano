"""ISR reweighting of the chargino signal (SIG-07) with its up/down variation (SYST-05), Run 2 method.

Legacy OSUT3Analysis ISRWeightProducer: the weight is the Z->mumu data/MC scale factor evaluated at the
generator-level pT of the electroweakino pair -- the vector sum of every chargino (1000024) and neutralino
(1000022) GenPart that is a last copy but not a first copy (requireLastNotFirstCopy, Pythia8 style; this
drops the neutralino from the chargino decay, which is first and last copy at once). Beyond the last bin the
last bin applies. Up/down = SF +- its statistical error, in every bin together.

Run 2 multiplied a second, MadGraph/PYTHIA factor because its signal was PYTHIA-only; the Run 3 signal is
madgraph-pythia8, so only the data/MC factor is used (decided 2026-10-15).

Scale factors come from limits/.../scripts/make_isr_scale_factors.py, one per MC campaign, picked by the
dataset's ``year`` metadata; the 2023D-template runs for 2024-2026 are 2023_postBPix and so use 2023D's
(decided 2026-10-23). The file is read once at configuration time and embedded, so Dask workers need no
extra transferred file.
"""

from __future__ import annotations

import json
from pathlib import Path

import awkward as ak
import numpy as np

from pocket_coffea.lib.weights.weights import WeightLambda

DEFAULT_PATH = Path(__file__).resolve().parent / "data" / "isr" / "isr_scale_factors_v1.json"
ERA_FOR_YEAR = {
    "2022_preEE": "2022CD",
    "2022_postEE": "2022EFG",
    "2023_preBPix": "2023C",
    "2023_postBPix": "2023D",
}
SUSY_PDGIDS = (1000022, 1000024)


def load_scale_factors(path=DEFAULT_PATH):
    table = json.loads(Path(path).read_text())
    return {era: table[era] for era in ERA_FOR_YEAR.values() if era in table}


def electroweakino_pair_pt(genpart):
    """Generator-level pT of the chargino/neutralino pair (ISRWeightProducer's isrPt)."""
    flags = genpart.statusFlags
    last = ((flags >> 13) & 1) == 1
    first = ((flags >> 12) & 1) == 1
    pid = abs(genpart.pdgId)
    susy = (pid == SUSY_PDGIDS[0]) | (pid == SUSY_PDGIDS[1])
    sel = susy & last & ~first
    px = ak.sum(ak.where(sel, genpart.pt * np.cos(genpart.phi), 0.0), axis=1)
    py = ak.sum(ak.where(sel, genpart.pt * np.sin(genpart.phi), 0.0), axis=1)
    return np.hypot(ak.to_numpy(px), ak.to_numpy(py))


def evaluate(scale_factors, pt):
    """(nominal, up, down) per event; edges[-1] is +inf, so pT beyond the last edge uses the last bin."""
    edges = np.asarray(scale_factors["edges"], dtype=float)
    sf = np.asarray(scale_factors["sf"], dtype=float)
    err = np.asarray(scale_factors["err"], dtype=float)
    idx = np.clip(np.searchsorted(edges, pt, side="right") - 1, 0, len(sf) - 1)
    return sf[idx], sf[idx] + err[idx], np.maximum(sf[idx] - err[idx], 0.0)


def make_isr_weight(scale_factors):
    """Return the ``isr`` weight class (has_variations: isrUp/isrDown)."""

    def compute(params, metadata, events, size, shape_variations):
        era = ERA_FOR_YEAR.get(metadata["year"])
        if era not in scale_factors:
            raise ValueError(f"isr weight: no scale factor for year={metadata['year']!r} (era {era!r})")
        return evaluate(scale_factors[era], electroweakino_pair_pt(events.GenPart))

    return WeightLambda.wrap_func(name="isr", function=compute, has_variations=True)


def make_leading_muon_trigger_sf_weight(collection="ISRZMuMuLeadingMuon"):
    """IsoMu24 trigger scale factor on the leading selected muon only (``sf_mu_trigger_leading``).

    PocketCoffea's ``sf_mu_trigger`` multiplies the per-muon SF over every muon in the collection, which is
    not right for a single-muon trigger in a two-muon selection; here only the leading Z muon is used, with
    PocketCoffea's own lookup (MUO-POG muon_Z.json, pT clipped to the 26 GeV map minimum). For the Z->mumu
    ISR control sample, as Run 2's weightsWithMuonSF applied a muon trigger SF there.
    """
    from pocket_coffea.lib.scale_factors import get_mu_sf

    def compute(params, metadata, events, size, shape_variations):
        muons = events[collection]
        pt = np.clip(ak.flatten(muons.pt), 26.0, 9999.0)
        sf, up, down = get_mu_sf(params, metadata["year"], pt, ak.flatten(muons.eta), ak.num(muons.pt), "trigger")
        return ak.prod(sf, axis=1), ak.prod(up, axis=1), ak.prod(down, axis=1)

    return WeightLambda.wrap_func(name="sf_mu_trigger_leading", function=compute, has_variations=True)
