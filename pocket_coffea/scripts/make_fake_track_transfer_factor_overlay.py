#!/usr/bin/env python3
"""Overlay several periods' dxy transfer-factor fits on one CMS-style plot.

Reuses the same folded-|dxy| fit (fit_dxy_transfer_factor) and signed-dxy
histogram that make-standard-fake-track-estimate's --fit-plots draws per
period, but layers every requested period onto a single axes -- e.g. a
Figure-7.15/Table-34-style comparison across run periods, styled like a
CMS-PAS overlay plot (CMS Preliminary + sqrt(s) label, ROOT-style ticks,
one color per period with its integrated luminosity in the legend).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_DIR = REPO_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from disapptrks.fake_tracks import fit_dxy_transfer_factor, plot_dxy_transfer_factor_overlay
from disapptrks.cli import _load_outputs, summed_hist_counts_edges  # noqa: E402  (script, not a package import)

CONTROL_CFG = {
    "zmumu": {
        "control_region": r"$Z\to\mu\mu$",
        "histogram": "fakeZMuMuFitTrack_absDxy",
        "signed_histogram": "fakeZMuMuFitTrack_dxy",
        "sample": "DATA_Muon",
    },
    "zee": {
        "control_region": r"$Z\to ee$",
        "histogram": "fakeZeeFitTrack_absDxy",
        "signed_histogram": "fakeZeeFitTrack_dxy",
        "sample": "DATA_EGamma",
    },
}

# pocket_coffea/parameters/lumi.yaml (PocketCoffea upstream), picobarns -> fb^-1,
# mapped from that file's Run-3 era names to this analysis's period names
# (2022CD == 2022_preEE, 2022EFG == 2022_postEE, 2023C == 2023_preBPix,
# 2023D == 2023_postBPix). 2026 has no lumi.yaml entry yet; the value below is
# the same user-provided full-Run-3-projection figure already used in
# disapptrks-datacards-limits (see that skill's own note on where it came from).
PERIOD_LUMI_FB = {
    "2022CD": 7980.4 / 1000.0,
    "2022EFG": 26671.7 / 1000.0,
    "2023C": 17960.0 / 1000.0,
    "2023D": 9680.0 / 1000.0,
    "2024": 109950.0 / 1000.0,
    "2025": 110840.0 / 1000.0,
    "2026": 25.31,
}

DEFAULT_PERIODS = ("2022CD", "2022EFG", "2023C", "2023D", "2024", "2025", "2026")


def _period_entry(period: str, coffea_file: Path, control: str) -> dict:
    cfg = CONTROL_CFG[control]
    outputs = _load_outputs([coffea_file])
    counts, edges = summed_hist_counts_edges(outputs, cfg["histogram"], sample=cfg["sample"])
    signed_counts, signed_edges = summed_hist_counts_edges(
        outputs, cfg["signed_histogram"], sample=cfg["sample"]
    )
    fit = fit_dxy_transfer_factor(
        counts,
        edges,
        control_region=cfg["control_region"],
        histogram=cfg["histogram"],
    )
    lumi = PERIOD_LUMI_FB.get(period)
    label = f"{period} ({lumi:.3g}" + r" fb$^{-1}$)" if lumi is not None else period
    return {
        "label": label,
        "counts": counts,
        "edges": edges,
        "signed_counts": signed_counts,
        "signed_edges": signed_edges,
        "fit": fit,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--control",
        choices=("zmumu", "zee"),
        required=True,
        help="Z control region to overlay",
    )
    parser.add_argument(
        "--periods",
        nargs="+",
        default=list(DEFAULT_PERIODS),
        help=f"Run periods to overlay (default: {' '.join(DEFAULT_PERIODS)})",
    )
    parser.add_argument(
        "--input-base",
        type=Path,
        default=Path("pocket_coffea/analysis_output"),
        help=(
            "Base containing <period>/fake_tracks/<control>_dedx/output_all.coffea "
            "(default: pocket_coffea/analysis_output)"
        ),
    )
    parser.add_argument(
        "--dedx-suffix",
        default="_dedx",
        help="Suffix on the per-control output directory (default: _dedx)",
    )
    parser.add_argument("-o", "--output", type=Path, required=True, help="Output plot path")
    parser.add_argument("--cms-label", default="CMS Preliminary")
    parser.add_argument("--energy-text", default="13.6 TeV")
    args = parser.parse_args()

    entries = []
    for period in args.periods:
        coffea_file = (
            args.input_base
            / period
            / "fake_tracks"
            / f"{args.control}{args.dedx_suffix}"
            / "output_all.coffea"
        )
        if not coffea_file.exists():
            raise SystemExit(f"error: missing {coffea_file} for period {period}")
        entries.append(_period_entry(period, coffea_file, args.control))

    plot_dxy_transfer_factor_overlay(
        entries,
        args.output,
        control_region=CONTROL_CFG[args.control]["control_region"],
        cms_label=args.cms_label,
        energy_text=args.energy_text,
    )
    print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
