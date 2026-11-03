"""Data vs hit-drop-corrected background MC (dissertation Fig. 7.23 style).

Applies each era's fitted probabilities (fit_results_v1.json from
fit_missing_hits_corrections.py) to the MC histograms, process by process,
with the same correction code the fit uses: missing middle hits with P0;
missing outer hits with the pre-TOB model (P_pre-TOB, P_TOB) for 0-5 and the
post-TOB model (P_post-TOB, P_TOB) for 6+, matching the two fit ranges as the
legacy findMissingOuterHitsCorrection.py assembled its corrected histogram.
MC is normalized to the data total; the uncorrected MC total is overlaid
(dashed) for reference.

Run from pocket_coffea/ inside ./shell:
    python scripts/plot_missing_hits_corrected.py [--eras "2022 CD" ...]
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
from coffea.util import load

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fit_missing_hits_corrections import (  # noqa: E402
    FITS,
    MIDDLE,
    NMIDDLE,
    NOUTER,
    OUTER,
    corrected_middle,
    corrected_outer,
)
from plot_missing_hits_data_mc import COLORS, GROUPS, LUMI_FB, plot  # noqa: E402

RESULTS = Path("analysis_output/missing_hits_corrections/fit_results_v1.json")
OUTDIR = Path("plots/missing_hits/corrected")


def _group_hist(output, name, samples):
    values = variances = None
    for sample in samples:
        by_dataset = output["variables"][name].get(sample, {})
        for hist in by_dataset.values():
            h = hist[{"cat": hist.axes["cat"][0]}]
            if "variation" in [a.name for a in h.axes]:
                h = h[{"variation": "nominal"}]
            v, w2 = h.values(), h.variances()
            values = v if values is None else values + v
            variances = w2 if variances is None else variances + w2
    return values, variances


def _data_1d(output, name):
    total = None
    for by_sample in output["variables"][name].values():
        for hist in by_sample.values():
            h = hist[{"cat": hist.axes["cat"][0]}]
            v = h.project(h.axes[0].name).values()
            total = v if total is None else total + v
    return total


def corrected_groups(mc_out, fit):
    """Per-group corrected (unnormalized) middle and outer distributions."""
    p0 = fit["middle"]["P0"][0]
    p_tob = fit["outer"]["P_TOB"][0]
    p_pre = fit["outer"]["P_preTOB"][0]
    p_post = fit["outer"]["P_postTOB"][0]
    middle, outer, middle_raw, outer_raw = [], [], [], []
    for _, samples in GROUPS:
        mid_v, mid_w2 = _group_hist(mc_out, MIDDLE, samples)
        out_v, out_w2 = _group_hist(mc_out, OUTER, samples)
        if mid_v is None:
            middle.append((None, None))
            outer.append((None, None))
            middle_raw.append(None)
            outer_raw.append(None)
            continue
        m, m_err = corrected_middle(p0, mid_v, mid_w2, normalize=False)
        pre, pre_err = corrected_outer(p_pre, p_tob, out_v, out_w2, "pre", normalize=False)
        post, post_err = corrected_outer(p_post, p_tob, out_v, out_w2, "post", normalize=False)
        o = np.where(np.arange(NOUTER) < 6, pre, post)
        o_err = np.where(np.arange(NOUTER) < 6, pre_err, post_err)
        middle.append((m, m_err**2))
        outer.append((o, o_err**2))
        middle_raw.append(mid_v.sum(axis=1)[:NMIDDLE])
        outer_raw.append(out_v.sum(axis=(1, 2))[:NOUTER])
    return middle, outer, middle_raw, outer_raw


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--eras", nargs="+", default=list(FITS))
    args = parser.parse_args()
    fits = json.loads(RESULTS.read_text())
    hep.style.use("CMS")
    OUTDIR.mkdir(parents=True, exist_ok=True)

    for era in args.eras:
        if era not in fits:
            print(f"--- {era}: no fit result, skipped")
            continue
        data_period, data_tag, mc_period, mc_tag = FITS[era]
        base = "analysis_output/{}/missing_hits_control/{}/output_all.coffea"
        data_out = load(base.format(data_period, data_tag))
        mc_out = load(base.format(mc_period, mc_tag))
        middle, outer, middle_raw, outer_raw = corrected_groups(mc_out, fits[era])

        fig, axes = plt.subplots(2, 2, figsize=(22, 11), sharex="col",
                                 gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05})
        print(f"=== {era}: data / corrected MC per bin (fractions)")
        for col, (name, groups, raw, nbins, xlabel) in enumerate((
            (MIDDLE, middle, middle_raw, NMIDDLE, "Missing middle hits"),
            (OUTER, outer, outer_raw, NOUTER, "Missing outer hits"),
        )):
            data = _data_1d(data_out, name)[:nbins]
            edges = np.arange(nbins + 1) - 0.5
            mc_total = sum(v.sum() for v, _ in groups if v is not None)
            stacked = plot(axes[0, col], axes[1, col], data, edges, groups,
                           data.sum() / mc_total, xlabel, nbins)
            raw_total = sum(r for r in raw if r is not None)
            axes[0, col].stairs(raw_total * data.sum() / raw_total.sum(), edges,
                                color="black", linestyle="--", linewidth=1.5,
                                label="Uncorrected MC")
            handles, labels = axes[0, col].get_legend_handles_labels()
            wanted = ["Data", *[g[0] for g in GROUPS], "MC stat. unc.", "Uncorrected MC"]
            order = [labels.index(w) for w in wanted if w in labels]
            axes[0, col].legend([handles[i] for i in order], [labels[i] for i in order],
                                fontsize=13, ncols=2, loc="upper right")
            ratios = [f"{data[i] / stacked[i]:.2f}" if stacked[i] > 0 else "-" for i in range(nbins)]
            print(f"   {xlabel}: " + " ".join(ratios))
            hep.cms.label("Preliminary", data=True, ax=axes[0, col], loc=0,
                          lumi=LUMI_FB.get(data_period), com=13.6, fontsize=18)
            axes[0, col].text(0.03, 0.95, f"{era} data vs\ncorrected {mc_period} MC\n(normalized to data)",
                              transform=axes[0, col].transAxes, fontsize=13, va="top")
            axes[1, col].set_ylim(0.01, 100)
        stem = f"missing_hits_corrected_{data_period}"
        for suffix in ("pdf", "png"):
            fig.savefig(OUTDIR / f"{stem}.{suffix}", bbox_inches="tight")
        plt.close(fig)
    print(f"Wrote plots to {OUTDIR}")


if __name__ == "__main__":
    main()
