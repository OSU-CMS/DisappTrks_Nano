from __future__ import annotations

import argparse
import json
from pathlib import Path

from . import greet
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

    args = parser.parse_args()
    if not hasattr(args, "func"):
        parser.print_help()
        return 0
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
