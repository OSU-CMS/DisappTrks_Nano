#!/usr/bin/env python3
"""Stage a subset of a PocketCoffea dataset JSON to OSU scratch.

This is intended for clusters where worker nodes can read
``/scratch0/user/$USER`` but cannot read the original long-term storage path.
It copies a slice of ROOT files into a scratch staging directory and writes a
new dataset JSON whose file paths point to the staged copies.
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path


DEFAULT_SCRATCH_BASE = Path("/scratch0/user/mjoyce/disapptrks_staging")


def flatten_dataset(dataset: dict) -> tuple[str, dict, list[str]]:
    if len(dataset) != 1:
        raise ValueError(
            "Expected exactly one dataset definition. "
            f"Found {len(dataset)}: {', '.join(dataset)}"
        )
    name, definition = next(iter(dataset.items()))
    files = list(definition.get("files", []))
    return name, definition, files


def safe_relative_source_path(path: Path) -> Path:
    """Return a relative path that preserves enough structure to avoid collisions."""

    parts = path.parts
    if "disapptrks_nano" in parts:
        index = parts.index("disapptrks_nano")
        return Path(*parts[index + 1 :])
    if "nano" in parts:
        index = parts.index("nano")
        return Path(*parts[index + 1 :])
    return Path(*parts[1:]) if path.is_absolute() else path


def copy_one(source: Path, destination: Path, *, force: bool, dry_run: bool) -> bool:
    if not source.is_file():
        raise FileNotFoundError(source)

    if destination.exists() and not force:
        if destination.stat().st_size == source.stat().st_size:
            return False
        raise FileExistsError(
            f"Destination exists with a different size: {destination}. "
            "Use --force to overwrite."
        )

    if dry_run:
        print(f"DRY-RUN copy {source} -> {destination}")
        return True

    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_suffix(destination.suffix + ".partial")
    if partial.exists():
        partial.unlink()
    shutil.copy2(source, partial)
    partial.replace(destination)
    return True


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Copy a slice of a dataset JSON to scratch and write a staged dataset JSON."
    )
    parser.add_argument("dataset_json", type=Path, help="Input dataset JSON.")
    parser.add_argument(
        "tag",
        help="Scratch subdirectory name, e.g. 2022CD_EGamma_electron_pveto_batch000.",
    )
    parser.add_argument(
        "-o",
        "--output-json",
        type=Path,
        help="Output staged dataset JSON. Default: datasets/staged_<tag>.json",
    )
    parser.add_argument(
        "--scratch-base",
        type=Path,
        default=DEFAULT_SCRATCH_BASE,
        help=f"Scratch staging parent. Default: {DEFAULT_SCRATCH_BASE}",
    )
    parser.add_argument(
        "--offset",
        type=int,
        default=0,
        help="First file index to stage.",
    )
    parser.add_argument(
        "--max-files",
        type=int,
        help="Maximum number of files to stage. Default: all files after offset.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite staged files if they already exist.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print planned copies and write no files.",
    )
    args = parser.parse_args()

    with args.dataset_json.open() as handle:
        dataset = json.load(handle)

    dataset_name, definition, files = flatten_dataset(dataset)
    if args.offset < 0:
        raise SystemExit("--offset must be non-negative")
    if args.offset >= len(files):
        raise SystemExit(f"--offset {args.offset} is outside 0..{len(files) - 1}")

    selected = files[args.offset :]
    if args.max_files is not None:
        if args.max_files <= 0:
            raise SystemExit("--max-files must be positive if set")
        selected = selected[: args.max_files]

    stage_dir = args.scratch_base / args.tag
    staged_files: list[str] = []
    n_copied = 0
    n_reused = 0
    for source_text in selected:
        source = Path(source_text)
        relative = safe_relative_source_path(source)
        destination = stage_dir / relative
        copied = copy_one(source, destination, force=args.force, dry_run=args.dry_run)
        if copied:
            n_copied += 1
        else:
            n_reused += 1
        staged_files.append(str(destination))

    output_json = args.output_json or Path("datasets") / f"staged_{args.tag}.json"
    staged_definition = dict(definition)
    staged_definition["files"] = staged_files
    staged_dataset = {dataset_name: staged_definition}

    if args.dry_run:
        print(f"DRY-RUN would write {output_json}")
    else:
        output_json.parent.mkdir(parents=True, exist_ok=True)
        with output_json.open("w") as handle:
            json.dump(staged_dataset, handle, indent=2, sort_keys=True)
            handle.write("\n")

    print(f"Input dataset:  {args.dataset_json}")
    print(f"Dataset name:   {dataset_name}")
    print(f"Selected files: {len(selected)} of {len(files)} starting at offset {args.offset}")
    print(f"Stage dir:      {stage_dir}")
    print(f"Output JSON:    {output_json}")
    print(f"Copied/reused:  {n_copied}/{n_reused}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
