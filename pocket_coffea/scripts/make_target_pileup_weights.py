"""Build pileup weights that reweight Summer23BPix MC to another year's data.

The LUM POG ``puWeights`` file for each year is data(year) / MC(that year's
campaign), so the official 2024/2025 weights assume a Summer24 MC pileup
profile and do not apply to Summer23BPix MC.  The 2024-2026 signal yields
(and the missing-hits correction fits for those years) use 2023D
(Summer23BPix) MC, which needs

    w(n) = P_data,target(n) / P_Summer23BPix(n).

The Summer23BPix profile is recovered exactly from the official 2023D
weights, which are data(2023D) / Summer23BPix:

    P_Summer23BPix(n) = P_data,2023D(n) / w_official,2023D(n)

(nominal 69.2 mb data profile).  Up/down use the target year's data profile
at 72.4/66.0 mb, the LUM convention.  2026 has no published data pileup
profile yet (2026-10-16), so its entry reuses the 2025 profile as a labelled
stand-in until LUM publishes one.

Inputs are in data/pileup/inputs/ (downloaded from
https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions2{3,4,5}/PileUp
and the LUM correctionlib payload on /cvmfs/cms-griddata.cern.ch).

Run from pocket_coffea/:
    python scripts/make_target_pileup_weights.py [--check-mc-filelist FILE]
"""

from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

import correctionlib
import correctionlib.schemav2 as cs
import numpy as np
import uproot

INPUTS = Path("data/pileup/inputs")
OUTPUT = Path("data/pileup/puWeights_Summer23BPix_to_target.json.gz")

OFFICIAL_2023D = INPUTS / "puWeights_Summer23BPix_2023D_official.json.gz"
OFFICIAL_2023D_NAME = "Collisions2023_369803_370790_eraD_GoldenJson"
DATA_2023D = (
    "pileupHistogram-Cert_Collisions2023_369803_370790_eraD_GoldenJson-"
    "13p6TeV-{xs}ub-99bins.root"
)
TARGET_DATA = {
    "2024": "dataPileupHistogram-2024BCDEFGHI_Golden-{xs}ub.root",
    "2025": "dataPileupHistogram-2025pp_Golden-{xs}ub.root",
}
# Target years without a published data profile, and the year standing in.
STAND_INS = {"2026": "2025"}
VARIATIONS = {"nominal": 69200, "up": 72400, "down": 66000}

NBINS = 100  # integer nTrueInt bins on [0, 100]
EDGES = np.arange(NBINS + 1, dtype=float)


def data_profile(filename: str) -> np.ndarray:
    """Data pileup histogram on the common [0, 100] integer binning."""
    hist = uproot.open(INPUTS / filename)["pileup"]
    values = hist.values()
    edges = hist.axis().edges()
    if not np.allclose(np.diff(edges), 1.0) or edges[0] != 0.0:
        raise ValueError(f"{filename}: expected unit bins starting at 0")
    profile = np.zeros(NBINS)
    count = min(len(values), NBINS)
    profile[:count] = values[:count]
    if len(values) > NBINS:
        profile[-1] += values[NBINS:].sum()
    return profile / profile.sum()


def summer23bpix_mc_profile() -> np.ndarray:
    weights = correctionlib.CorrectionSet.from_file(str(OFFICIAL_2023D))[
        OFFICIAL_2023D_NAME
    ]
    centers = EDGES[:-1] + 0.5
    official = np.array([weights.evaluate(n, "nominal") for n in centers])
    data = data_profile(DATA_2023D.format(xs=VARIATIONS["nominal"]))
    mc = np.divide(data, official, out=np.zeros(NBINS), where=official > 0)
    return mc / mc.sum()


def ratio(numerator: np.ndarray, denominator: np.ndarray) -> np.ndarray:
    return np.divide(
        numerator, denominator, out=np.zeros(NBINS), where=denominator > 0
    )


def target_correction(target: str, source: str, mc: np.ndarray) -> cs.Correction:
    if target != source:
        profile = f"STAND-IN: {source} profile, no {target} profile published yet"
    else:
        profile = "LUM golden-JSON profile"
    description = (
        f"Pileup weight reweighting Summer23BPix MC to the {target} data "
        f"pileup profile ({profile}); MC profile = P_data,2023D / official "
        "2023D puWeights; up/down = 72.4/66.0 mb min-bias cross section."
    )
    content = []
    for variation, xs in VARIATIONS.items():
        values = ratio(data_profile(TARGET_DATA[source].format(xs=xs)), mc)
        content.append(
            cs.CategoryItem(
                key=variation,
                value=cs.Binning(
                    nodetype="binning",
                    input="NumTrueInteractions",
                    edges=list(EDGES),
                    content=[float(v) for v in values],
                    flow="clamp",
                ),
            )
        )
    return cs.Correction(
        name=f"Summer23BPix_to_{target}",
        description=description,
        version=1,
        inputs=[
            cs.Variable(name="NumTrueInteractions", type="real", description="Pileup_nTrueInt"),
            cs.Variable(name="weights", type="string", description="nominal, up, or down"),
        ],
        output=cs.Variable(name="weight", type="real"),
        data=cs.Category(nodetype="category", input="weights", content=content),
    )


def check_mc_profile(mc: np.ndarray, filelist: Path) -> None:
    """Compare the derived MC profile with Pileup_nTrueInt in real MC files."""
    counts = np.zeros(NBINS)
    files = [line.strip() for line in filelist.read_text().splitlines() if line.strip()]
    for path in files:
        with uproot.open(path) as handle:
            events = handle["Events"]
            ntrue = events["Pileup_nTrueInt"].array(library="np")
            sign = np.sign(events["genWeight"].array(library="np"))
        counts += np.histogram(np.clip(ntrue, 0, NBINS - 1e-6), bins=EDGES, weights=sign)[0]
    observed = counts / counts.sum()
    mean_derived = float((EDGES[:-1] + 0.5) @ mc)
    mean_observed = float((EDGES[:-1] + 0.5) @ observed)
    print(f"MC check: {len(files)} file(s), {counts.sum():.0f} signed events")
    print(f"  mean nTrueInt derived {mean_derived:.2f}  observed {mean_observed:.2f}")
    populated = observed > 1e-3
    pulls = (observed - mc)[populated] / np.sqrt(observed[populated] / counts.sum())
    print(f"  bins with >0.1% of events: {populated.sum()}, max |pull| {np.abs(pulls).max():.1f}")
    for lo in range(0, NBINS, 10):
        print(f"  n={lo:2d}-{lo + 9:2d}  derived {mc[lo:lo + 10].sum():.4f}  observed {observed[lo:lo + 10].sum():.4f}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--check-mc-filelist", type=Path)
    args = parser.parse_args()

    mc = summer23bpix_mc_profile()
    print(f"Summer23BPix MC profile: mean nTrueInt {(EDGES[:-1] + 0.5) @ mc:.2f}, "
          f"empty bins {[i for i in range(NBINS) if mc[i] == 0]}")

    corrections = []
    for target in [*TARGET_DATA, *STAND_INS]:
        source = STAND_INS.get(target, target)
        correction = target_correction(target, source, mc)
        corrections.append(correction)
        nominal = data_profile(TARGET_DATA[source].format(xs=VARIATIONS["nominal"]))
        lost = nominal[mc == 0].sum()
        print(f"{correction.name}: data mean nTrueInt {(EDGES[:-1] + 0.5) @ nominal:.2f}, "
              f"data fraction in bins with no MC {lost:.2e}")

    cset = cs.CorrectionSet(schema_version=2, corrections=corrections)
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(OUTPUT, "wt") as handle:
        handle.write(cset.model_dump_json(exclude_unset=True) if hasattr(cset, "model_dump_json") else cset.json(exclude_unset=True))
    print(f"Wrote {OUTPUT}")

    # Closure: Summer23BPix MC weighted by each target weight must reproduce
    # that target's data profile.
    loaded = correctionlib.CorrectionSet.from_file(str(OUTPUT))
    centers = EDGES[:-1] + 0.5
    for target in [*TARGET_DATA, *STAND_INS]:
        source = STAND_INS.get(target, target)
        weights = np.array([loaded[f"Summer23BPix_to_{target}"].evaluate(n, "nominal") for n in centers])
        reweighted = mc * weights
        nominal = data_profile(TARGET_DATA[source].format(xs=VARIATIONS["nominal"]))
        print(f"closure {target}: max |reweighted MC - data| per bin "
              f"{np.abs(reweighted / reweighted.sum() - nominal).max():.2e}, "
              f"mean weight over MC {reweighted.sum():.4f}")

    if args.check_mc_filelist:
        check_mc_profile(mc, args.check_mc_filelist)


if __name__ == "__main__":
    main()
