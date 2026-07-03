#!/usr/bin/env python3
"""Copy ROOT files listed in OSUNano filelists to another cluster.

Each input text file creates one destination directory named after the filelist
stem.  For local destinations, files are copied directly into the destination.
For remote destinations, files are first copied from XRootD to a local staging
area with ``xrdcp`` and then synchronized to the destination with ``rsync``.
The EOS path below a configurable prefix is preserved to avoid collisions
between many files named ``nano_1.root``.

Example
-------

  python3 scripts/copy_filelists_to_remote.py \
    filelists/Muon_2023C.txt filelists/EGamma_2023C.txt \
    --dest mjoyce@cms-t3.mps.ohio-state.edu:/abyss/users/mjoyce/disapptrks_nano

The spelling ``user@host://path`` is accepted and normalized to the rsync/scp
form ``user@host:/path``.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


DEFAULT_STRIP_PREFIX = "/store/group/lpcdisapptrks/nano"


def read_filelist(path: Path) -> list[str]:
    files = []
    for line in path.read_text().splitlines():
        item = line.strip()
        if not item or item.startswith("#"):
            continue
        files.append(item)
    return files


def normalize_remote_base(dest: str) -> str:
    """Normalize a remote destination to rsync/scp's ``host:/path`` form."""
    if "://" in dest and not dest.startswith(("root://", "http://", "https://")):
        host, path = dest.split("://", 1)
        return f"{host}:/{path.lstrip('/')}"
    return dest.rstrip("/")


def split_remote(dest: str) -> tuple[str, str] | None:
    """Return ``(host, path)`` for rsync-style remote destinations."""
    if ":" not in dest:
        return None
    host, path = dest.split(":", 1)
    if "/" not in path:
        return None
    return host, path.rstrip("/")


def url_path(url: str) -> str:
    """Extract the file path from an XRootD URL or return a local-looking path."""
    if "://" in url:
        scheme_sep = url.find("://")
        first_slash = url.find("/", scheme_sep + 3)
        path = url[first_slash:] if first_slash >= 0 else ""
    else:
        path = url

    while path.startswith("//"):
        path = path[1:]
    if not path.startswith("/"):
        path = "/" + path
    return path


def relative_output_path(url: str, strip_prefix: str) -> Path:
    path = url_path(url)
    prefix = strip_prefix.rstrip("/")
    if path == prefix:
        rel = Path(path).name
    elif path.startswith(prefix + "/"):
        rel = path[len(prefix) + 1 :]
    else:
        rel = path.lstrip("/")
    return Path(rel)


def run(cmd: list[str], *, dry_run: bool = False) -> None:
    print("+ " + " ".join(cmd), flush=True)
    if not dry_run:
        subprocess.run(cmd, check=True)


def copy_one(
    source: str,
    destination: Path,
    *,
    dry_run: bool = False,
    overwrite: bool = False,
) -> Path:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and not overwrite:
        print(f"skip existing {destination}", flush=True)
        return destination

    tmp_destination = destination.with_suffix(destination.suffix + ".partial")
    if tmp_destination.exists():
        tmp_destination.unlink()

    run(["xrdcp", "-f", source, str(tmp_destination)], dry_run=dry_run)
    if not dry_run:
        tmp_destination.replace(destination)
    return destination


def sync_stage(stage_dir: Path, dest: str, *, dry_run: bool = False) -> None:
    remote = split_remote(dest)
    if remote is None:
        Path(dest).mkdir(parents=True, exist_ok=True)
        run(["rsync", "-a", f"{stage_dir}/", f"{dest}/"], dry_run=dry_run)
        return

    host, remote_path = remote
    run(["ssh", host, "mkdir", "-p", remote_path], dry_run=dry_run)
    run(["rsync", "-a", f"{stage_dir}/", f"{dest}/"], dry_run=dry_run)


def process_filelist(
    filelist: Path,
    *,
    dest_base: str,
    stage_base: Path,
    strip_prefix: str,
    max_files: int | None,
    workers: int,
    dry_run: bool,
    overwrite: bool,
    keep_stage: bool,
    continue_on_error: bool,
) -> int:
    sources = read_filelist(filelist)
    if max_files is not None:
        sources = sources[:max_files]

    name = filelist.stem
    is_remote_destination = split_remote(dest_base) is not None
    stage_dir = stage_base / name if is_remote_destination else Path(dest_base) / name
    if is_remote_destination and stage_dir.exists() and overwrite:
        shutil.rmtree(stage_dir)
    stage_dir.mkdir(parents=True, exist_ok=True)

    mode = "stage+rsync" if is_remote_destination else "direct local copy"
    print(
        f"\n== {filelist} -> {dest_base}/{name} ({len(sources)} file(s), {mode}) ==",
        flush=True,
    )

    failures = 0

    def _copy(source: str) -> Path:
        rel = relative_output_path(source, strip_prefix)
        return copy_one(
            source,
            stage_dir / rel,
            dry_run=dry_run,
            overwrite=overwrite,
        )

    if workers <= 1:
        for source in sources:
            try:
                _copy(source)
            except Exception as exc:  # noqa: BLE001 - continue-on-error reporting
                failures += 1
                print(f"ERROR copying {source}: {exc}", file=sys.stderr, flush=True)
                if not continue_on_error:
                    raise
    else:
        with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
            future_to_source = {pool.submit(_copy, source): source for source in sources}
            for future in concurrent.futures.as_completed(future_to_source):
                source = future_to_source[future]
                try:
                    future.result()
                except Exception as exc:  # noqa: BLE001 - continue-on-error reporting
                    failures += 1
                    print(f"ERROR copying {source}: {exc}", file=sys.stderr, flush=True)
                    if not continue_on_error:
                        raise

    if is_remote_destination:
        sync_stage(stage_dir, f"{dest_base}/{name}", dry_run=dry_run)

    if is_remote_destination and not keep_stage and not dry_run:
        shutil.rmtree(stage_dir)

    return failures


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Copy ROOT files from one or more filelists to a remote cluster."
    )
    parser.add_argument("filelists", nargs="+", type=Path)
    parser.add_argument(
        "--dest",
        required=True,
        help="Destination base, e.g. user@host:/abyss/users/me/disapptrks_nano.",
    )
    parser.add_argument(
        "--stage-dir",
        type=Path,
        help="Local staging directory. Default: create a temporary directory.",
    )
    parser.add_argument(
        "--strip-prefix",
        default=DEFAULT_STRIP_PREFIX,
        help=(
            "Prefix stripped from EOS paths before preserving subdirectories. "
            f"Default: {DEFAULT_STRIP_PREFIX}"
        ),
    )
    parser.add_argument("--max-files", type=int, help="Only copy the first N files per list.")
    parser.add_argument("--workers", type=int, default=4, help="Parallel xrdcp workers.")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite staged files.")
    parser.add_argument("--keep-stage", action="store_true", help="Do not delete staged files.")
    parser.add_argument(
        "--continue-on-error",
        action="store_true",
        help="Continue copying other files if one xrdcp fails.",
    )
    parser.add_argument("--dry-run", action="store_true", help="Print commands only.")
    args = parser.parse_args(argv)

    dest_base = normalize_remote_base(args.dest)
    missing = [str(path) for path in args.filelists if not path.exists()]
    if missing:
        raise SystemExit("Missing filelist(s): " + ", ".join(missing))

    if args.stage_dir is None:
        tmp = tempfile.TemporaryDirectory(prefix="disapptrks_nano_copy_")
        stage_base = Path(tmp.name)
    else:
        tmp = None
        stage_base = args.stage_dir
        stage_base.mkdir(parents=True, exist_ok=True)

    try:
        failures = 0
        for filelist in args.filelists:
            failures += process_filelist(
                filelist,
                dest_base=dest_base,
                stage_base=stage_base,
                strip_prefix=args.strip_prefix,
                max_files=args.max_files,
                workers=args.workers,
                dry_run=args.dry_run,
                overwrite=args.overwrite,
                keep_stage=args.keep_stage,
                continue_on_error=args.continue_on_error,
            )
    finally:
        if tmp is not None:
            tmp.cleanup()

    if failures:
        print(f"Completed with {failures} failed copy operation(s).", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
