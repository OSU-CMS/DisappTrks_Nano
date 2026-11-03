"""Write per-campaign background-MC dataset JSONs for the custom NanoAOD.

One JSON per Run 3 MC campaign (eos_<period>_BkgMC_OSUv3.json), holding the
legacy Run 2-style background mix used by the missing-hits correction fits
(W->lnu, Z->ll, ttbar, QCD, single top, diboson, Z->nunu). Each physics
process is one dataset entry: an ``_ext1`` extension is merged with its
nominal sample (PocketCoffea normalizes every dataset entry to xsec*lumi by
its own sum of generator weights, so separate entries would double count).
``sample`` groups processes as in the legacy plot legend.

Cross sections (pb) are the legacy Run 3 values from
OSUT3Analysis/Configuration/python/configurationOptions.py (crossSections),
the ones the dissertation used; they are identical across 2022-2023
campaigns there. Overlapping QCD samples (the flat QCD_PT-15to7000 and the
Summer23BPix-only QCD_PT-600) are skipped.

``nevents`` is the exact Events-entry count (PocketCoffea normalizes MC by
the sum of generator weights, not by this, but it is the reference for
checking that a run processed every event).

Run from pocket_coffea/ inside ./shell (needs xrdfs and uproot):
    python scripts/make_bkg_mc_dataset_jsons.py
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
from pathlib import Path

from disapptrks.datasets import count_root_events

XROOTD = "root://cmseos.fnal.gov"
BASE = "/store/group/lpcdisapptrks/nano/dev_v2"
NANO_TAG = "OSUv3"

# period label -> (MC campaign, dataset-metadata year, era)
CAMPAIGNS = {
    "2022CD": ("Summer22", "2022_preEE", "CD"),
    "2022EFG": ("Summer22EE", "2022_postEE", "EFG"),
    "2023C": ("Summer23", "2023_preBPix", "C"),
    "2023D": ("Summer23BPix", "2023_postBPix", "D"),
}

TUNE = "_TuneCP5_13p6TeV"
# primary dataset -> (sample group, xsec in pb)
PROCESSES = {
    f"WtoLNu-4Jets{TUNE}_madgraphMLM-pythia8": ("MC_WtoLNu", 64451.4),
    f"DYto2L-4Jets_MLL-50{TUNE}_madgraphMLM-pythia8": ("MC_DYto2L", 6331.5),
    f"DYto2L-4Jets_MLL-10to50{TUNE}_madgraphMLM-pythia8": ("MC_DYto2L", 20950.0),
    f"TTto2L2Nu{TUNE}_powheg-pythia8": ("MC_TT", 98.88),
    f"TTtoLNu2Q{TUNE}_powheg-pythia8": ("MC_TT", 407.43),
    f"TTto4Q{TUNE}_powheg-pythia8": ("MC_TT", 419.69),
    f"TBbartoLplusNuBbar-s-channel-4FS{TUNE}_amcatnlo-pythia8": ("MC_SingleTop", 2.30),
    f"TbarBtoLminusNuB-s-channel-4FS{TUNE}_amcatnlo-pythia8": ("MC_SingleTop", 1.35),
    f"TQbartoLNu-t-channel{TUNE}_powheg-pythia8": ("MC_SingleTop", 43.91),
    f"TbarQtoLNu-t-channel{TUNE}_powheg-pythia8": ("MC_SingleTop", 26.18),
    f"TWminusto2L2Nu{TUNE}_powheg-pythia8": ("MC_SingleTop", 14.38),
    f"TWminustoLNu2Q{TUNE}_powheg-pythia8": ("MC_SingleTop", 29.63),
    f"TbarWplusto2L2Nu{TUNE}_powheg-pythia8": ("MC_SingleTop", 14.38),
    f"TbarWplustoLNu2Q{TUNE}_powheg-pythia8": ("MC_SingleTop", 29.63),
    f"WW{TUNE}_pythia8": ("MC_Diboson", 80.23),
    f"WZ{TUNE}_pythia8": ("MC_Diboson", 29.1),
    f"ZZ{TUNE}_pythia8": ("MC_Diboson", 12.75),
    f"Zto2Nu-4Jets_HT-100to200{TUNE}_madgraphMLM-pythia8": ("MC_Zto2Nu", 314.03),
    f"Zto2Nu-4Jets_HT-200to400{TUNE}_madgraphMLM-pythia8": ("MC_Zto2Nu", 87.15),
    f"Zto2Nu-4Jets_HT-400to800{TUNE}_madgraphMLM-pythia8": ("MC_Zto2Nu", 15.13),
    f"Zto2Nu-4Jets_HT-800to1500{TUNE}_madgraphMLM-pythia8": ("MC_Zto2Nu", 1.56),
    f"Zto2Nu-4Jets_HT-1500to2500{TUNE}_madgraphMLM-pythia8": ("MC_Zto2Nu", 0.113),
    f"Zto2Nu-4Jets_HT-2500{TUNE}_madgraphMLM-pythia8": ("MC_Zto2Nu", 0.0077),
    f"QCD_PT-15to30{TUNE}_pythia8": ("MC_QCD", 1301000000.0),
    f"QCD_PT-30to50{TUNE}_pythia8": ("MC_QCD", 113300000.0),
    f"QCD_PT-50to80{TUNE}_pythia8": ("MC_QCD", 16760000.0),
    f"QCD_PT-80to120{TUNE}_pythia8": ("MC_QCD", 2534000.0),
    f"QCD_PT-120to170{TUNE}_pythia8": ("MC_QCD", 445800.0),
    f"QCD_PT-170to300{TUNE}_pythia8": ("MC_QCD", 113700.0),
    f"QCD_PT-300to470{TUNE}_pythia8": ("MC_QCD", 7559.0),
    f"QCD_PT-470to600{TUNE}_pythia8": ("MC_QCD", 626.4),
    f"QCD_PT-600to800{TUNE}_pythia8": ("MC_QCD", 178.6),
    f"QCD_PT-800to1000{TUNE}_pythia8": ("MC_QCD", 30.57),
    f"QCD_PT-1000to1400{TUNE}_pythia8": ("MC_QCD", 8.92),
    f"QCD_PT-1400to1800{TUNE}_pythia8": ("MC_QCD", 0.8103),
    f"QCD_PT-1800to2400{TUNE}_pythia8": ("MC_QCD", 0.1148),
    f"QCD_PT-2400to3200{TUNE}_pythia8": ("MC_QCD", 0.007542),
    f"QCD_PT-3200{TUNE}_pythia8": ("MC_QCD", 0.0002331),
}


def xrdfs_root_files(path: str) -> list[str]:
    result = subprocess.run(
        ["xrdfs", XROOTD, "ls", "-R", "-u", path],
        check=True,
        capture_output=True,
        text=True,
    )
    return sorted(
        line.strip()
        for line in result.stdout.splitlines()
        if line.strip().endswith(".root") and "/failed/" not in line
    )


def task_dirs(primary: str, campaign: str) -> list[str]:
    """Nominal and extension production directories for one campaign."""
    listing = subprocess.run(
        ["xrdfs", XROOTD, "ls", f"{BASE}/{primary}"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.split()
    pattern = re.compile(
        rf"_Run3{campaign}MiniAODv4(_ext\d+)?_customNanoAOD_{NANO_TAG}$"
    )
    return [path for path in listing if pattern.search(path)]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--periods", nargs="+", default=list(CAMPAIGNS))
    parser.add_argument("--event-count-workers", type=int, default=24)
    args = parser.parse_args()

    for period in args.periods:
        campaign, year, era = CAMPAIGNS[period]
        datasets = {}
        for primary, (sample, xsec) in PROCESSES.items():
            dirs = task_dirs(primary, campaign)
            files = [f for d in dirs for f in xrdfs_root_files(d)]
            if not files:
                print(f"  WARNING {period}: no files for {primary}")
                continue
            nevents, _ = count_root_events(files, max_workers=args.event_count_workers)
            name = f"{primary.replace(TUNE, '')}_{campaign}_{NANO_TAG}"
            datasets[name] = {
                "files": files,
                "metadata": {
                    "sample": sample,
                    "year": year,
                    "era": era,
                    "isMC": "True",
                    "primaryDataset": primary,
                    "nevents": str(nevents),
                    "nano_version": 12,
                    "xsec": str(xsec),
                },
            }
            print(f"  {period} {name}: {len(files)} files from {len(dirs)} "
                  f"task dir(s), {nevents} events", flush=True)
        output = Path("datasets") / f"eos_{period}_BkgMC_{NANO_TAG}.json"
        output.write_text(json.dumps(datasets, indent=2) + "\n")
        print(f"Wrote {output}: {len(datasets)} datasets, "
              f"{sum(len(d['files']) for d in datasets.values())} files")


if __name__ == "__main__":
    main()
