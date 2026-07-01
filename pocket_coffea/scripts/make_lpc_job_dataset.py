#!/usr/bin/env python3
"""Build a per-job PocketCoffea dataset JSON from a larger dataset JSON.

This helper is intentionally dependency-free so it can run inside the LPC
Condor sandbox before the PocketCoffea job starts.  It preserves the original
dataset metadata and keeps only the file slice assigned to the Condor ProcId.
"""

from __future__ import annotations

import argparse
import copy
import json
import math
from pathlib import Path


def count_files(dataset: dict) -> int:
    return sum(len(definition.get("files", [])) for definition in dataset.values())


def slice_dataset(
    dataset: dict,
    job_id: int,
    files_per_job: int,
    fallback_events_per_file: int,
) -> dict:
    start = job_id * files_per_job
    stop = start + files_per_job

    output = {}
    cursor = 0
    for name, definition in dataset.items():
        files = list(definition.get("files", []))
        next_cursor = cursor + len(files)

        local_start = max(start - cursor, 0)
        local_stop = min(stop - cursor, len(files))
        selected = files[local_start:local_stop] if local_start < local_stop else []

        if selected:
            new_definition = copy.deepcopy(definition)
            new_definition["files"] = selected
            metadata = dict(new_definition.get("metadata", {}))
            try:
                nevents = int(metadata.get("nevents", "0"))
            except (TypeError, ValueError):
                nevents = 0
            if nevents <= 0:
                nevents = max(len(selected) * fallback_events_per_file, 1)
            metadata["nevents"] = str(nevents)
            new_definition["metadata"] = metadata
            output[name] = new_definition

        cursor = next_cursor

    return output


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--job-id", required=True, type=int)
    parser.add_argument("--files-per-job", required=True, type=int)
    parser.add_argument("--fallback-events-per-file", type=int, default=50000)
    args = parser.parse_args()

    if args.files_per_job <= 0:
        raise SystemExit("--files-per-job must be positive")

    dataset = json.loads(args.input.read_text())
    n_files = count_files(dataset)
    n_jobs = math.ceil(n_files / args.files_per_job) if n_files else 0

    selected = slice_dataset(
        dataset,
        args.job_id,
        args.files_per_job,
        args.fallback_events_per_file,
    )
    if not selected:
        raise SystemExit(
            f"job {args.job_id} has no files; dataset has {n_files} files "
            f"and {n_jobs} jobs for files_per_job={args.files_per_job}"
        )

    args.output.write_text(json.dumps(selected, indent=2) + "\n")
    print(
        f"Wrote {args.output} with {count_files(selected)} file(s) "
        f"for job {args.job_id}/{n_jobs}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
