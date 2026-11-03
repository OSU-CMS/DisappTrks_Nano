"""Cumulative missing-hit fractions: data vs uncorrected vs hit-drop-corrected MC.

For each period: fraction of muon-control-region tracks with >= n missing outer
hits (n = 1..4) and fraction of tracks (no inner/middle requirement) with 0
missing middle hits, in data, uncorrected MC and corrected MC (fit_results_v1.json
probabilities, corrected exactly as in plot_missing_hits_corrected.py). These are
the quantities the selection (>= 3 missing outer, 0 missing middle) and the Run 2
missing-hits systematics depend on.

Run from pocket_coffea/ inside ./shell:
    python scripts/missing_hits_cumulative.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from coffea.util import load

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fit_missing_hits_corrections import (  # noqa: E402
    FITS,
    MIDDLE,
    NMIDDLE,
    NOUTER,
    OUTER,
    _summed,
    corrected_middle,
    corrected_outer,
)

RESULTS = Path("analysis_output/missing_hits_corrections/fit_results_v1.json")
OUTPUT = Path("analysis_output/missing_hits_corrections/cumulative_fractions_v1.json")
THRESHOLDS = (1, 2, 3, 4)


def fraction(values, variances, select):
    """Fraction in ``select`` with its uncertainty (binomial-like for weighted counts)."""
    total = values.sum()
    f = values[select].sum() / total
    # Var(f) for weighted counts: sum over bins of (df/dv_i)^2 var_i
    deriv = np.where(select, (1 - f) / total, -f / total)
    return f, float(np.sqrt((deriv**2 * variances).sum()))


def main():
    fits = json.loads(RESULTS.read_text())
    base = "analysis_output/{}/missing_hits_control/{}/output_all.coffea"
    out = {}
    for era, (data_period, data_tag, mc_period, mc_tag) in FITS.items():
        if era not in fits:
            continue
        fit = fits[era]
        data = load(base.format(data_period, data_tag))
        mc = load(base.format(mc_period, mc_tag))

        d_out, _ = _summed(data, OUTER)
        d_out = d_out.sum(axis=(1, 2))[:NOUTER]
        d_mid, _ = _summed(data, MIDDLE)
        d_mid = d_mid.sum(axis=1)
        mc_out_v, mc_out_w2 = _summed(mc, OUTER)
        mc_mid_v, mc_mid_w2 = _summed(mc, MIDDLE)

        raw_out = (mc_out_v.sum(axis=(1, 2))[:NOUTER], mc_out_w2.sum(axis=(1, 2))[:NOUTER])
        p_tob = fit["outer"]["P_TOB"][0]
        pre, pre_err = corrected_outer(fit["outer"]["P_preTOB"][0], p_tob, mc_out_v, mc_out_w2,
                                       "pre", normalize=False)
        post, post_err = corrected_outer(fit["outer"]["P_postTOB"][0], p_tob, mc_out_v, mc_out_w2,
                                         "post", normalize=False)
        low = np.arange(NOUTER) < 6
        cor_out = (np.where(low, pre, post), np.where(low, pre_err, post_err) ** 2)

        raw_mid = (mc_mid_v.sum(axis=1)[:NMIDDLE], mc_mid_w2.sum(axis=1)[:NMIDDLE])
        cm, cm_err = corrected_middle(fit["middle"]["P0"][0], mc_mid_v, mc_mid_w2, normalize=False)
        cor_mid = (cm, cm_err**2)

        bins_out = np.arange(NOUTER)
        bins_mid = np.arange(NMIDDLE)
        rows = {}
        for n in THRESHOLDS:
            sel = bins_out >= n
            rows[f"outer>={n}"] = {
                "data": fraction(d_out, d_out, sel),
                "mc_uncorrected": fraction(*raw_out, sel),
                "mc_corrected": fraction(*cor_out, sel),
            }
        rows["middle==0"] = {
            "data": fraction(d_mid[:NMIDDLE], d_mid[:NMIDDLE], bins_mid == 0),
            "mc_uncorrected": fraction(*raw_mid, bins_mid == 0),
            "mc_corrected": fraction(*cor_mid, bins_mid == 0),
        }
        out[era] = rows

        print(f"=== {era}")
        print(f"   {'quantity':11s} {'data':>20s} {'MC uncorrected':>20s} {'MC corrected':>20s}"
              f" {'data/unc':>9s} {'data/cor':>9s}")
        for key, r in rows.items():
            (fd, ed), (fu, eu), (fc, ec) = r["data"], r["mc_uncorrected"], r["mc_corrected"]
            print(f"   {key:11s} {fd:.3e} ± {ed:.1e}   {fu:.3e} ± {eu:.1e}   {fc:.3e} ± {ec:.1e}"
                  f" {fd / fu:9.3f} {fd / fc:9.3f}")

    OUTPUT.write_text(json.dumps(out, indent=2) + "\n")
    print(f"Wrote {OUTPUT}")


if __name__ == "__main__":
    main()
