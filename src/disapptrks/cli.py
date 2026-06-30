from __future__ import annotations

import argparse
import json
from pathlib import Path

from . import greet
from .datasets import (
    build_dataset_definition,
    list_eos_root_files,
    root_files_from_lines,
    write_dataset_definition,
)
from .schema import audit_root_file
from .summaries import (
    summarize_ss_subtracted_veto_probability,
    summarize_veto_probability,
)


def _audit_command(args: argparse.Namespace) -> int:
    report = audit_root_file(args.file)
    print(json.dumps(report.as_dict(), indent=2, sort_keys=True))
    return 0 if report.ready_for(args.scope) else 2


def _summarize_pveto_command(args: argparse.Namespace) -> int:
    from coffea.util import load

    output = load(args.file)
    if args.ss_subtract:
        summary = summarize_ss_subtracted_veto_probability(
            output["cutflow"],
            os_denominator_name=args.denominator,
            os_numerator_name=args.numerator,
            ss_denominator_name=args.ss_denominator,
            ss_numerator_name=args.ss_numerator,
            dataset=args.dataset,
            sample=args.sample,
            variation=args.variation,
        )

        if args.json:
            print(json.dumps(summary.as_dict(), indent=2, sort_keys=True))
            return 0

        print(f"OS denominator ({summary.os_denominator_name}): {summary.os_denominator:g}")
        print(f"OS numerator   ({summary.os_numerator_name}): {summary.os_numerator:g}")
        print(f"SS denominator ({summary.ss_denominator_name}): {summary.ss_denominator:g}")
        print(f"SS numerator   ({summary.ss_numerator_name}): {summary.ss_numerator:g}")
        print(f"Subtracted denominator: {summary.denominator:g}")
        print(f"Subtracted numerator:   {summary.numerator:g}")
        print(f"P(pass veto): {summary.probability:.6g} ± {summary.uncertainty:.6g}")
        return 0

    summary = summarize_veto_probability(
        output["cutflow"],
        denominator_name=args.denominator,
        numerator_name=args.numerator,
        dataset=args.dataset,
        sample=args.sample,
        variation=args.variation,
    )

    if args.json:
        print(json.dumps(summary.as_dict(), indent=2, sort_keys=True))
        return 0

    print(f"Denominator ({summary.denominator_name}): {summary.denominator:g}")
    print(f"Numerator   ({summary.numerator_name}): {summary.numerator:g}")
    print(f"P(pass veto): {summary.probability:.6g} ± {summary.uncertainty:.6g}")
    return 0


def _make_dataset_json_command(args: argparse.Namespace) -> int:
    if args.filelist:
        files = root_files_from_lines(args.filelist.read_text().splitlines())
    else:
        if not args.eos_path:
            raise SystemExit("error: eos_path is required unless --filelist is used")
        files = list_eos_root_files(
            args.eos_path,
            xrootd=args.xrootd,
            recursive=args.recursive,
        )

    dataset = build_dataset_definition(
        dataset_name=args.dataset_name,
        files=files,
        sample=args.sample,
        year=args.year,
        era=args.era,
        primary_dataset=args.primary_dataset,
        is_mc=args.is_mc,
        nevents=args.nevents,
        nano_version=args.nano_version,
    )
    write_dataset_definition(dataset, args.output)

    print(f"Wrote {len(files)} ROOT files to {args.output}")
    if not files:
        print("Warning: no ROOT files found")
        return 2
    return 0


def main():
    parser = argparse.ArgumentParser(prog="disapptrks")
    subparsers = parser.add_subparsers(dest="command")

    audit = subparsers.add_parser(
        "audit-schema",
        help="Check a custom NanoAOD file against the analysis branch contract.",
    )
    audit.add_argument("file", type=Path)
    audit.add_argument(
        "--scope",
        choices=("search", "backgrounds", "fiducial-maps"),
        default="backgrounds",
    )
    audit.set_defaults(func=_audit_command)

    pveto = subparsers.add_parser(
        "summarize-pveto",
        help="Summarize a veto pass probability from a PocketCoffea output file.",
    )
    pveto.add_argument("file", type=Path)
    pveto.add_argument("--dataset", help="Restrict to one dataset key.")
    pveto.add_argument("--sample", help="Restrict to one sample key.")
    pveto.add_argument("--variation", default="nominal")
    pveto.add_argument("--json", action="store_true", help="Print JSON output.")
    pveto.add_argument(
        "--ss-subtract",
        action="store_true",
        help="Compute (OS numerator - SS numerator) / (OS denominator - SS denominator).",
    )
    pveto.add_argument(
        "--denominator",
        default="muon_veto_zwindow",
        help="OS cutflow category used as denominator.",
    )
    pveto.add_argument(
        "--numerator",
        default="muon_veto_zwindow_pass",
        help="OS cutflow category used as numerator.",
    )
    pveto.add_argument(
        "--ss-denominator",
        default="muon_veto_ss_zwindow",
        help="SS cutflow category subtracted from the denominator.",
    )
    pveto.add_argument(
        "--ss-numerator",
        default="muon_veto_ss_zwindow_pass",
        help="SS cutflow category subtracted from the numerator.",
    )
    pveto.set_defaults(func=_summarize_pveto_command)

    dataset_json = subparsers.add_parser(
        "make-dataset-json",
        help="Create a PocketCoffea dataset JSON from an EOS directory or filelist.",
    )
    dataset_json.add_argument("eos_path", nargs="?", help="EOS directory, e.g. /store/user/...")
    dataset_json.add_argument("-o", "--output", type=Path, required=True)
    dataset_json.add_argument("--dataset-name", required=True)
    dataset_json.add_argument("--sample", default="DATA_Muon")
    dataset_json.add_argument("--year", required=True)
    dataset_json.add_argument("--era", required=True)
    dataset_json.add_argument("--primary-dataset", default="Muon")
    dataset_json.add_argument("--nevents", default="0")
    dataset_json.add_argument("--nano-version", type=int, default=15)
    dataset_json.add_argument("--is-mc", action="store_true")
    dataset_json.add_argument(
        "--xrootd",
        default="root://cmseos.fnal.gov",
        help="XRootD endpoint passed to xrdfs.",
    )
    dataset_json.add_argument(
        "-R",
        "--recursive",
        action="store_true",
        help="Use xrdfs ls -R -u to recurse through subdirectories.",
    )
    dataset_json.add_argument(
        "--filelist",
        type=Path,
        help="Read ROOT file URLs from a text file instead of calling xrdfs.",
    )
    dataset_json.set_defaults(func=_make_dataset_json_command)

    args = parser.parse_args()
    if not hasattr(args, "func"):
        parser.print_help()
        return 0
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
