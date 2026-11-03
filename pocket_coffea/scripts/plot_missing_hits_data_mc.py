"""Uncorrected data vs background-MC missing-hit distributions (Fig. 7.22 style).

For each comparison: data in the muon control region vs the stacked,
uncorrected background MC (no hit-drop correction), for missing middle hits
(control region without the missing inner/middle requirements) and missing
outer hits (full muon control region). MC is normalized to the data total of
each distribution; the raw MC/data normalization is printed. 2024 and 2025
compare to 2023D MC reweighted to that year's pileup profile.

Run from pocket_coffea/ inside ./shell:
    python scripts/plot_missing_hits_data_mc.py
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
from coffea.util import load

OUTDIR = Path("plots/missing_hits/uncorrected")
MIDDLE = "missingHitsNoInnerMiddle_missingMiddle_stripLayers"
OUTER = "missingHitsMuonCtrl_missingOuter_TOBLayers_stripLayers"

# data period, MC period, MC output tag, label
COMPARISONS = [
    ("2022CD", "2022CD", "bkgmc_v1", "2022CD data vs Summer22 MC"),
    ("2022EFG", "2022EFG", "bkgmc_v1", "2022EFG data vs Summer22EE MC"),
    ("2023C", "2023C", "bkgmc_v1", "2023C data vs Summer23 MC"),
    ("2023D", "2023D", "bkgmc_v1", "2023D data vs Summer23BPix MC"),
    ("2024", "2023D", "bkgmc_v1_pu2024", "2024 data vs Summer23BPix MC (2024 PU)"),
    ("2025", "2023D", "bkgmc_v1_pu2025", "2025 data vs Summer23BPix MC (2025 PU)"),
    ("2026", "2023D", "bkgmc_v1_pu2026", "2026 data vs Summer23BPix MC (2025 PU stand-in)"),
]
LUMI_FB = {"2022CD": 7.98, "2022EFG": 26.67, "2023C": 18.08, "2023D": 9.69}

# Legend order as in dissertation Fig. 7.22/D.1; stacked in reverse (smallest
# at the bottom, W on top) so small processes stay visible on the log axis.
# Fixed colors from the CMS 10-color palette.
GROUPS = [
    ("W$\\to\\ell\\nu$", ("MC_WtoLNu",)),
    ("Z$\\to\\ell\\ell$", ("MC_DYto2L",)),
    ("$t\\bar{t}$", ("MC_TT",)),
    ("QCD", ("MC_QCD",)),
    ("Single top", ("MC_SingleTop",)),
    ("Diboson", ("MC_Diboson",)),
    ("Z$\\to\\nu\\bar{\\nu}$", ("MC_Zto2Nu",)),
]
COLORS = ["#3f90da", "#ffa90e", "#bd1f01", "#832db6", "#94a4a2", "#e76300", "#717581"]
# Dissertation Fig. D.1 y-axis; eras with ~10x more events get a higher top.
Y_RANGE = (1e-3, 5e6)
Y_RANGE_HIGH_STATS = (1e-3, 5e7)


def _projected(hist, axis_index):
    hist = hist[{"cat": hist.axes["cat"][0]}]
    if "variation" in [axis.name for axis in hist.axes]:
        hist = hist[{"variation": "nominal"}]
    projected = hist.project(hist.axes[axis_index].name)
    return projected.values(), projected.variances(), projected.axes[0].edges


def _leaves(node):
    if isinstance(node, dict):
        for child in node.values():
            yield from _leaves(child)
    else:
        yield node


def data_distribution(output, name):
    values = edges = None
    for hist in _leaves(output["variables"][name]):
        v, _, e = _projected(hist, 0)
        values = v if values is None else values + v
        edges = e
    return values, edges


def mc_distributions(output, name):
    by_group = []
    for _, samples in GROUPS:
        values = variances = None
        for sample in samples:
            for hist in _leaves(output["variables"][name].get(sample, {})):
                v, w2, _ = _projected(hist, 0)
                values = v if values is None else values + v
                variances = w2 if variances is None else variances + w2
        by_group.append((values, variances))
    return by_group


def plot(ax_main, ax_ratio, data, edges, groups, scale, xlabel, nbins):
    data, edges = data[:nbins], edges[: nbins + 1]
    stacked = np.zeros(nbins)
    total_var = np.zeros(nbins)
    for (label, _), color, (values, variances) in reversed(list(zip(GROUPS, COLORS, groups))):
        if values is None:
            continue
        v, w2 = values[:nbins] * scale, variances[:nbins] * scale**2
        ax_main.stairs(stacked + v, edges, baseline=stacked, fill=True, color=color,
                       label=label, edgecolor="white", linewidth=0.6)
        stacked += v
        total_var += w2
    mc_err = np.sqrt(total_var)
    ax_main.stairs(stacked + mc_err, edges, baseline=stacked - mc_err, fill=True,
                   facecolor="none", edgecolor="black", hatch="////", linewidth=0,
                   label="MC stat. unc.")
    centers = 0.5 * (edges[1:] + edges[:-1])
    nonzero = data > 0
    ax_main.errorbar(centers[nonzero], data[nonzero], yerr=np.sqrt(data[nonzero]),
                     fmt="o", color="black", markersize=5, label="Data")
    ax_main.set_yscale("log")
    ax_main.set_ylim(*(Y_RANGE if data.sum() < 1e6 else Y_RANGE_HIGH_STATS))
    ax_main.set_ylabel("Events")
    handles, labels = ax_main.get_legend_handles_labels()
    order = [labels.index("Data")] + [labels.index(g[0]) for g in GROUPS if g[0] in labels] + [labels.index("MC stat. unc.")]
    ax_main.legend([handles[i] for i in order], [labels[i] for i in order],
                   fontsize=13, ncols=2, loc="upper right")

    ratio_mask = stacked > 0
    ratio = np.full(nbins, np.nan)
    ratio[ratio_mask] = data[ratio_mask] / stacked[ratio_mask]
    ratio_err = np.full(nbins, np.nan)
    ratio_err[ratio_mask] = np.sqrt(data[ratio_mask]) / stacked[ratio_mask]
    rel_mc = np.where(ratio_mask, mc_err / np.where(ratio_mask, stacked, 1), 0)
    ax_ratio.stairs(1 + rel_mc, edges, baseline=1 - rel_mc, fill=True, facecolor="none",
                    edgecolor="black", hatch="////", linewidth=0)
    ax_ratio.errorbar(centers, ratio, yerr=ratio_err, fmt="o", color="black", markersize=5)
    ax_ratio.axhline(1, color="0.4", linewidth=1)
    ax_ratio.set_yscale("log")
    ax_ratio.set_ylim(0.3, 300)
    ax_ratio.set_ylabel("Data / MC")
    ax_ratio.set_xlabel(xlabel)
    ax_ratio.set_xlim(edges[0], edges[-1])
    return stacked


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--periods", nargs="+", help="data periods to plot (default: all)")
    parser.add_argument("--data-tag", default="data_v1")
    parser.add_argument("--mc-suffix", default="",
                        help="appended to each MC output tag, e.g. _noHPdEdx")
    parser.add_argument("--label-suffix", default="", help="appended to the plot label")
    parser.add_argument("--file-suffix", default="", help="appended to the output file name")
    args = parser.parse_args()

    hep.style.use("CMS")
    OUTDIR.mkdir(parents=True, exist_ok=True)
    for data_period, mc_period, tag, label in COMPARISONS:
        if args.periods and data_period not in args.periods:
            continue
        tag += args.mc_suffix
        label += args.label_suffix
        data_out = load(f"analysis_output/{data_period}/missing_hits_control/{args.data_tag}/output_all.coffea")
        mc_out = load(f"analysis_output/{mc_period}/missing_hits_control/{tag}/output_all.coffea")
        fig, axes = plt.subplots(2, 2, figsize=(22, 11), sharex="col",
                                 gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05})
        print(f"=== {label}")
        for col, (name, xlabel, nbins) in enumerate(
            ((MIDDLE, "Missing middle hits", 6), (OUTER, "Missing outer hits", 16))
        ):
            data, edges = data_distribution(data_out, name)
            groups = mc_distributions(mc_out, name)
            mc_total = sum(v.sum() for v, _ in groups if v is not None)
            scale = data.sum() / mc_total
            if data_period == mc_period:
                print(f"   {xlabel}: raw MC/data normalization {1 / scale:.3f}")
            stacked = plot(axes[0, col], axes[1, col], data, edges, groups, scale, xlabel, nbins)
            data_frac = data[:nbins] / data.sum()
            mc_frac = stacked / stacked.sum()
            print(f"   {xlabel} fraction per bin, data / MC:")
            for i in range(nbins):
                if data_frac[i] > 0 or mc_frac[i] > 0:
                    print(f"      {i:2d}: {data_frac[i]:.3e} / {mc_frac[i]:.3e}"
                          f"  ratio {data_frac[i] / mc_frac[i] if mc_frac[i] else float('nan'):.2f}")
            hep.cms.label("Preliminary", data=True, ax=axes[0, col], loc=0,
                          lumi=LUMI_FB.get(data_period), com=13.6, fontsize=18)
            axes[0, col].text(0.03, 0.95, label + "\nuncorrected MC, normalized to data", va="top",
                              transform=axes[0, col].transAxes, fontsize=13)
        for suffix in ("pdf", "png"):
            fig.savefig(OUTDIR / f"missing_hits_uncorrected_{data_period}{args.file_suffix}.{suffix}", bbox_inches="tight")
        plt.close(fig)
    print(f"Wrote plots to {OUTDIR}")


if __name__ == "__main__":
    main()
