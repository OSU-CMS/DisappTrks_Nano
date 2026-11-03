"""Missing inner / middle hit systematics (SYST-06/07), published Run 2 method.

EXO-19-010 / legacy HitsSystematic driver ``HITS``: in the muon control region
without the missing inner/middle requirements (missing_hits_control mode,
category missing_hits_ctrl_no_inner_middle), compare the fraction of tracks
with 0 missing middle hits (eff_X) and with 0 missing inner hits (eff_Y)
between data and simulation:

    systematic = |eff_data - eff_MC| / eff_data,

one value per period, applied to every layer bin. Middle hits use the
hit-drop-corrected MC (fit_results_v1.json probabilities, analytic correction
validated against the random correction); inner hits are not corrected, so
they use the MC as is. 2024-2026 use the 2023D MC reweighted to that year's
pileup, as in the fits. Muon control region only (no HITS_V2), per Matt
2026-10-22.

Run from pocket_coffea/ inside ./shell:
    python scripts/compute_hits_systematics.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from coffea.util import load

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fit_missing_hits_corrections import FITS, MIDDLE, NMIDDLE, _summed, corrected_middle  # noqa: E402

MIDDLE_INNER = "missingHitsNoInnerMiddle_missingMiddle_missingInner"  # [missing middle, missing inner]
RESULTS = Path("analysis_output/missing_hits_corrections/fit_results_v1.json")
OUTPUT = Path("analysis_output/missing_hits_corrections/hits_systematics_v1.json")
PERIOD = {"2022 CD": "2022CD", "2022 EFG": "2022EFG", "2023 C": "2023C", "2023 D": "2023D",
          "2024": "2024", "2025": "2025", "2026": "2026"}


def zero_fraction(values, variances):
    """Fraction in bin 0 and its uncertainty for (weighted) counts."""
    total = values.sum()
    f = values[0] / total
    deriv = np.full(len(values), -f / total)
    deriv[0] = (1 - f) / total
    return float(f), float(np.sqrt((deriv**2 * variances).sum()))


def systematic(data, mc):
    (fd, ed), (fm, em) = data, mc
    value = abs(fd - fm) / fd
    error = np.hypot(em / fd, fm * ed / fd**2)  # propagated statistical uncertainty
    return 100 * value, 100 * error


def main():
    fits = json.loads(RESULTS.read_text())
    base = "analysis_output/{}/missing_hits_control/{}/output_all.coffea"
    out = {}
    print(f"{'period':8s} {'eff0 middle data':>18s} {'MC corr.':>18s} {'syst %':>12s}"
          f"   {'eff0 inner data':>18s} {'MC':>18s} {'syst %':>12s}")
    for era, (data_period, data_tag, mc_period, mc_tag) in FITS.items():
        data = load(base.format(data_period, data_tag))
        mc = load(base.format(mc_period, mc_tag))
        d2, _ = _summed(data, MIDDLE_INNER)
        m2, m2_w2 = _summed(mc, MIDDLE_INNER)
        mid_v, mid_w2 = _summed(mc, MIDDLE)

        d_middle = d2.sum(axis=1)
        d_inner = d2.sum(axis=0)
        corr, corr_err = corrected_middle(fits[era]["middle"]["P0"][0], mid_v, mid_w2, normalize=False)
        middle_data = zero_fraction(d_middle, d_middle)
        middle_mc = zero_fraction(corr, corr_err**2)
        inner_data = zero_fraction(d_inner, d_inner)
        inner_mc = zero_fraction(m2.sum(axis=0), m2_w2.sum(axis=0))
        s_mid, s_mid_err = systematic(middle_data, middle_mc)
        s_in, s_in_err = systematic(inner_data, inner_mc)
        period = PERIOD[era]
        out[period] = {
            "missing_middle_hits": {"percent": s_mid, "stat_err_percent": s_mid_err,
                                    "eff_data": middle_data, "eff_mc_corrected": middle_mc},
            "missing_inner_hits": {"percent": s_in, "stat_err_percent": s_in_err,
                                   "eff_data": inner_data, "eff_mc": inner_mc},
            "inputs": {"data": f"{data_period}/{data_tag}", "mc": f"{mc_period}/{mc_tag}"},
        }
        print(f"{period:8s} {middle_data[0]:.5f}±{middle_data[1]:.5f} {middle_mc[0]:.5f}±{middle_mc[1]:.5f}"
              f" {s_mid:5.2f}±{s_mid_err:4.2f}   {inner_data[0]:.5f}±{inner_data[1]:.5f}"
              f" {inner_mc[0]:.5f}±{inner_mc[1]:.5f} {s_in:5.2f}±{s_in_err:4.2f}")
    OUTPUT.write_text(json.dumps(out, indent=2) + "\n")
    print(f"Wrote {OUTPUT}")


if __name__ == "__main__":
    main()
