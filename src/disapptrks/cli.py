from __future__ import annotations

import argparse
import json
from pathlib import Path

from . import greet
from .datasets import (
    build_dataset_definition,
    group_osunano_files,
    list_eos_root_files,
    root_files_from_lines,
    scan_eos_bases_for_root_files,
    write_dataset_definition,
    write_grouped_filelists,
)
from .fake_tracks import estimate_fake_track_background, write_fake_track_latex
from .schema import audit_root_file
from .summaries import (
    summarize_ss_subtracted_veto_probability,
    summarize_veto_probability,
)
from .tables import write_muon_cutflow_latex, write_muon_pveto_latex


def _sum_nested_numeric(left, right):
    if isinstance(left, dict) and isinstance(right, dict):
        out = dict(left)
        for key, value in right.items():
            out[key] = _sum_nested_numeric(out[key], value) if key in out else value
        return out
    if isinstance(left, (int, float)) and isinstance(right, (int, float)):
        return left + right
    return right


def _load_merged_cutflow(files: list[Path]) -> dict:
    from coffea.util import load

    merged = {}
    for path in files:
        output = load(path)
        merged = _sum_nested_numeric(merged, output["cutflow"])
    return merged


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


def _estimate_fake_tracks_command(args: argparse.Namespace) -> int:
    if args.counts_json:
        source = json.loads(args.counts_json.read_text())
        source_is_cutflow = False
    else:
        if not args.files:
            raise SystemExit("error: at least one coffea file is required unless --counts-json is used")
        source = _load_merged_cutflow(args.files)
        source_is_cutflow = True

    estimates = [
        estimate_fake_track_background(
            source,
            layer=layer,
            transfer_signal_category=args.transfer_signal_category,
            transfer_sideband_category=args.transfer_sideband_category,
            control_category=args.control_category,
            basic_yield_category=args.basic_yield_category,
            z_to_ll_yield_category=args.z_to_ll_yield_category,
            dataset=args.dataset,
            sample=args.sample,
            variation=args.variation,
            source_is_cutflow=source_is_cutflow,
            prescale=args.prescale,
        )
        for layer in args.layers
    ]

    if args.output_json:
        args.output_json.parent.mkdir(parents=True, exist_ok=True)
        args.output_json.write_text(json.dumps([e.as_dict() for e in estimates], indent=2, sort_keys=True))
        print(f"Wrote {args.output_json}")

    if args.output_tex:
        write_fake_track_latex(
            estimates,
            args.output_tex,
            run_period=args.run_period,
            include_table_env=args.table_env,
        )
        print(f"Wrote {args.output_tex}")

    for estimate in estimates:
        print(
            f"{estimate.layer}: "
            f"xi={estimate.transfer_factor.value:.6g} ± {estimate.transfer_factor.error:.6g}, "
            f"N_ctrl={estimate.control.value:.6g} ± {estimate.control.error:.6g}, "
            f"N_fake={estimate.estimate.value:.6g} ± {estimate.estimate.error:.6g}"
        )
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


def _make_era_filelists_command(args: argparse.Namespace) -> int:
    if args.filelist:
        files = root_files_from_lines(args.filelist.read_text().splitlines())
    else:
        files = scan_eos_bases_for_root_files(args.eos_paths, xrootd=args.xrootd)

    grouped = group_osunano_files(
        files,
        prod_version_policy=args.prod_version_policy,
    )
    outputs = write_grouped_filelists(
        grouped,
        output_dir=args.output_dir,
        dataset_json_dir=args.dataset_json_dir,
        nano_version=args.nano_version,
    )

    print(f"Scanned {len(files)} ROOT file(s)")
    for key, entry in outputs.items():
        suffix = ""
        if "dataset_json" in entry:
            suffix = f", json={entry['dataset_json']}"
        print(f"{key:18s} {entry['n_files']:6d} files -> {entry['filelist']}{suffix}")

    empty = [key for key, entry in outputs.items() if entry["n_files"] == 0]
    if empty:
        print("Warning: empty groups: " + ", ".join(empty))
        return 2 if args.fail_on_empty else 0
    return 0


def _make_pveto_tables_command(args: argparse.Namespace) -> int:
    from coffea.util import load

    output = load(args.file)
    cutflow = output["cutflow"]

    write_muon_cutflow_latex(
        cutflow,
        args.cutflow_tex,
        dataset=args.dataset,
        sample=args.sample,
        variation=args.variation,
        include_table_env=args.table_env,
    )
    summaries = write_muon_pveto_latex(
        cutflow,
        args.pveto_tex,
        run_period=args.run_period,
        flavor=args.flavor,
        layers=args.layers,
        dataset=args.dataset,
        sample=args.sample,
        variation=args.variation,
        os_denominator_name=args.denominator,
        os_numerator_name=args.numerator,
        ss_denominator_name=args.ss_denominator,
        ss_numerator_name=args.ss_numerator,
        include_table_env=args.table_env,
    )

    print(f"Wrote {args.cutflow_tex}")
    print(f"Wrote {args.pveto_tex}")
    for layer, summary in summaries.items():
        print(
            f"{layer}: P_veto = "
            f"{summary.central:.6g} +{summary.err_up:.6g} -{summary.err_down:.6g} "
            f"(signed numerator={summary.numerator:g}, denominator={summary.denominator:g})"
        )
    return 0


def _make_lepton_pveto_table_command(args: argparse.Namespace) -> int:
    cutflow = _load_merged_cutflow(args.files)

    if args.mode == "tau_mu":
        defaults = {
            "os_denominator_name": "tau_mu_veto_masswindow",
            "os_numerator_name": "tau_mu_pveto_masswindow_pass",
            "ss_denominator_name": "tau_mu_veto_ss_masswindow",
            "ss_numerator_name": "tau_mu_pveto_ss_masswindow_pass",
            "flavor": r"$\tau_{\mu}$",
        }
    elif args.mode == "tau_ele":
        defaults = {
            "os_denominator_name": "tau_ele_veto_masswindow",
            "os_numerator_name": "tau_ele_pveto_masswindow_pass",
            "ss_denominator_name": "tau_ele_veto_ss_masswindow",
            "ss_numerator_name": "tau_ele_pveto_ss_masswindow_pass",
            "flavor": r"$\tau_{e}$",
        }
    elif args.mode == "electron":
        defaults = {
            "os_denominator_name": "electron_veto_zwindow",
            "os_numerator_name": "electron_pveto_zwindow_pass",
            "ss_denominator_name": "electron_veto_ss_zwindow",
            "ss_numerator_name": "electron_pveto_ss_zwindow_pass",
            "flavor": r"$e$",
        }
    else:
        raise ValueError(f"unknown lepton Pveto mode: {args.mode}")

    summaries = write_muon_pveto_latex(
        cutflow,
        args.output,
        run_period=args.run_period,
        flavor=args.flavor or defaults["flavor"],
        layers=args.layers,
        dataset=args.dataset,
        sample=args.sample,
        variation=args.variation,
        os_denominator_name=args.denominator or defaults["os_denominator_name"],
        os_numerator_name=args.numerator or defaults["os_numerator_name"],
        ss_denominator_name=args.ss_denominator or defaults["ss_denominator_name"],
        ss_numerator_name=args.ss_numerator or defaults["ss_numerator_name"],
        include_table_env=args.table_env,
    )

    print(f"Wrote {args.output}")
    for layer, summary in summaries.items():
        print(
            f"{layer}: P_veto = "
            f"{summary.central:.6g} +{summary.err_up:.6g} -{summary.err_down:.6g} "
            f"(signed numerator={summary.numerator:g}, denominator={summary.denominator:g})"
        )
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
        default="muon_pveto_zwindow_pass",
        help="OS cutflow category used as numerator.",
    )
    pveto.add_argument(
        "--ss-denominator",
        default="muon_veto_ss_zwindow",
        help="SS cutflow category subtracted from the denominator.",
    )
    pveto.add_argument(
        "--ss-numerator",
        default="muon_pveto_ss_zwindow_pass",
        help="SS cutflow category subtracted from the numerator.",
    )
    pveto.set_defaults(func=_summarize_pveto_command)

    pveto_tables = subparsers.add_parser(
        "make-pveto-tables",
        help="Write AN-style muon Pveto cutflow and probability LaTeX tables.",
    )
    pveto_tables.add_argument("file", type=Path)
    pveto_tables.add_argument("--dataset", help="Restrict to one dataset key.")
    pveto_tables.add_argument("--sample", help="Restrict to one sample key.")
    pveto_tables.add_argument("--variation", default="nominal")
    pveto_tables.add_argument("--run-period", required=True, help="Run-period label used in the Pveto table.")
    pveto_tables.add_argument("--flavor", default=r"$\mu$", help="Flavor label used in the Pveto table.")
    pveto_tables.add_argument(
        "--layers",
        nargs="+",
        default=["NLayers4", "NLayers5", "NLayers6plus", "combinedBins"],
        help="Layer-bin rows to include in the Pveto table.",
    )
    pveto_tables.add_argument(
        "--cutflow-tex",
        type=Path,
        default=Path("muon_pveto_cutflow.tex"),
        help="Output path for the cutflow LaTeX table.",
    )
    pveto_tables.add_argument(
        "--pveto-tex",
        type=Path,
        default=Path("muon_pveto_table.tex"),
        help="Output path for the Pveto LaTeX table.",
    )
    pveto_tables.add_argument(
        "--table-env",
        action="store_true",
        help="Wrap each tabular in a LaTeX table environment.",
    )
    pveto_tables.add_argument(
        "--denominator",
        default="muon_veto_zwindow",
        help="OS cutflow category used as denominator.",
    )
    pveto_tables.add_argument(
        "--numerator",
        default="muon_pveto_zwindow_pass",
        help="OS cutflow category used as numerator.",
    )
    pveto_tables.add_argument(
        "--ss-denominator",
        default="muon_veto_ss_zwindow",
        help="SS cutflow category subtracted from the denominator.",
    )
    pveto_tables.add_argument(
        "--ss-numerator",
        default="muon_pveto_ss_zwindow_pass",
        help="SS cutflow category subtracted from the numerator.",
    )
    pveto_tables.set_defaults(func=_make_pveto_tables_command)

    lepton_pveto_table = subparsers.add_parser(
        "make-lepton-pveto-table",
        help=(
            "Write an AN-style Pveto probability table for electron or tau "
            "tag-and-probe outputs, merging one or more coffea files."
        ),
    )
    lepton_pveto_table.add_argument("files", nargs="+", type=Path)
    lepton_pveto_table.add_argument(
        "--mode",
        required=True,
        choices=("electron", "tau_mu", "tau_ele"),
        help="Category-name preset to use.",
    )
    lepton_pveto_table.add_argument("--dataset", help="Restrict to one dataset key.")
    lepton_pveto_table.add_argument("--sample", help="Restrict to one sample key.")
    lepton_pveto_table.add_argument("--variation", default="nominal")
    lepton_pveto_table.add_argument(
        "--run-period",
        required=True,
        help="Run-period label used in the Pveto table.",
    )
    lepton_pveto_table.add_argument(
        "--flavor",
        help="Override the flavor label used in the Pveto table.",
    )
    lepton_pveto_table.add_argument(
        "--layers",
        nargs="+",
        default=["NLayers4", "NLayers5", "NLayers6plus", "combinedBins"],
        help="Layer-bin rows to include in the Pveto table.",
    )
    lepton_pveto_table.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="Output path for the Pveto LaTeX table.",
    )
    lepton_pveto_table.add_argument(
        "--table-env",
        action="store_true",
        help="Wrap the tabular in a LaTeX table environment.",
    )
    lepton_pveto_table.add_argument(
        "--denominator",
        help="Override the OS cutflow category used as denominator.",
    )
    lepton_pveto_table.add_argument(
        "--numerator",
        help="Override the OS cutflow category used as numerator.",
    )
    lepton_pveto_table.add_argument(
        "--ss-denominator",
        help="Override the SS cutflow category subtracted from the denominator.",
    )
    lepton_pveto_table.add_argument(
        "--ss-numerator",
        help="Override the SS cutflow category subtracted from the numerator.",
    )
    lepton_pveto_table.set_defaults(func=_make_lepton_pveto_table_command)

    fake_tracks = subparsers.add_parser(
        "estimate-fake-tracks",
        help="Compute the AN-style fake-track estimate from coffea cutflow counts.",
    )
    fake_tracks.add_argument(
        "files",
        nargs="*",
        type=Path,
        help="PocketCoffea output files. Their cutflows are summed before estimating.",
    )
    fake_tracks.add_argument(
        "--counts-json",
        type=Path,
        help="Use a flat JSON mapping of category names to counts instead of coffea files.",
    )
    fake_tracks.add_argument("--dataset", help="Restrict to one dataset key for coffea input.")
    fake_tracks.add_argument("--sample", help="Restrict to one sample key for coffea input.")
    fake_tracks.add_argument("--variation", default="nominal")
    fake_tracks.add_argument("--run-period", required=True, help="Run-period label used in the LaTeX table.")
    fake_tracks.add_argument(
        "--layers",
        nargs="+",
        default=["NLayers4", "NLayers5", "NLayers6plus", "combinedBins"],
        help="Layer-bin rows to estimate. Category patterns may use {layer}.",
    )
    fake_tracks.add_argument(
        "--transfer-signal-category",
        default="fake_basic3hits_d0_signal",
        help="Category/count for the transfer-factor numerator, usually |d0| < 0.02 in the 3-hit/basic sample.",
    )
    fake_tracks.add_argument(
        "--transfer-sideband-category",
        default="fake_basic3hits_d0_sideband",
        help="Category/count for the transfer-factor denominator, usually the sideband 0.05 < |d0| < 0.5.",
    )
    fake_tracks.add_argument(
        "--control-category",
        default="fake_control_{layer}",
        help="Category/count for the target-layer inverted-d0 control yield. May use {layer}.",
    )
    fake_tracks.add_argument(
        "--basic-yield-category",
        help="Optional BasicSelection yield category for normalizing Z->ll to the search sample.",
    )
    fake_tracks.add_argument(
        "--z-to-ll-yield-category",
        help="Optional inclusive Z->ll yield category for normalizing Z->ll to the search sample.",
    )
    fake_tracks.add_argument(
        "--prescale",
        type=float,
        default=1.0,
        help="Optional prescale/completion factor applied to N_ctrl.",
    )
    fake_tracks.add_argument(
        "--output-json",
        type=Path,
        help="Write detailed estimate components to JSON.",
    )
    fake_tracks.add_argument(
        "--output-tex",
        type=Path,
        help="Write an AN-style LaTeX summary table.",
    )
    fake_tracks.add_argument(
        "--table-env",
        action="store_true",
        help="Wrap the LaTeX tabular in a table environment.",
    )
    fake_tracks.set_defaults(func=_estimate_fake_tracks_command)

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

    era_filelists = subparsers.add_parser(
        "make-era-filelists",
        help="Scan OSUNano EOS areas and write Muon/EGamma filelists split by Run-3 era group.",
    )
    era_filelists.add_argument(
        "eos_paths",
        nargs="*",
        default=[
            "/store/group/lpcdisapptrks/nano/dev",
            "/store/group/lpcdisapptrks/nano/prod",
        ],
        help="EOS base directories to scan recursively.",
    )
    era_filelists.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        default=Path("filelists"),
        help="Directory where .txt filelists are written.",
    )
    era_filelists.add_argument(
        "--dataset-json-dir",
        type=Path,
        help="Also write PocketCoffea dataset JSONs to this directory.",
    )
    era_filelists.add_argument(
        "--filelist",
        type=Path,
        help="Classify ROOT files from an existing text file instead of querying EOS.",
    )
    era_filelists.add_argument(
        "--xrootd",
        default="root://cmseos.fnal.gov",
        help="XRootD endpoint passed to xrdfs.",
    )
    era_filelists.add_argument(
        "--prod-version-policy",
        choices=("latest", "all"),
        default="all",
        help="How to handle prod directories with *_vN suffixes. Default treats versions as literal names and keeps all.",
    )
    era_filelists.add_argument(
        "--nano-version",
        type=int,
        default=None,
        help=(
            "Override the per-era NanoAOD version metadata. "
            "By default this uses 12 for 2022/2023 and 15 for 2024/2025."
        ),
    )
    era_filelists.add_argument(
        "--fail-on-empty",
        action="store_true",
        help="Return nonzero if any primary/year group has zero files.",
    )
    era_filelists.set_defaults(func=_make_era_filelists_command)

    args = parser.parse_args()
    if not hasattr(args, "func"):
        parser.print_help()
        return 0
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
