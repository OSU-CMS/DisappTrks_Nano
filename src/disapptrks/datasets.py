"""Helpers for building PocketCoffea dataset JSON files."""

from __future__ import annotations

import json
import subprocess
from pathlib import Path


def root_files_from_lines(lines: list[str]) -> list[str]:
    """Return sorted ROOT-file paths from command/filelist lines."""
    files = []
    for line in lines:
        path = line.strip()
        if not path or path.startswith("#"):
            continue
        if path.endswith(".root"):
            files.append(path)
    return sorted(files)


def list_eos_root_files(
    eos_path: str,
    *,
    xrootd: str = "root://cmseos.fnal.gov",
    recursive: bool = False,
) -> list[str]:
    """List ROOT files under an EOS path using ``xrdfs ls -u``."""
    command = ["xrdfs", xrootd, "ls", "-u"]
    if recursive:
        command.append("-R")
    command.append(eos_path)

    result = subprocess.run(
        command,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    return root_files_from_lines(result.stdout.splitlines())


def build_dataset_definition(
    *,
    dataset_name: str,
    files: list[str],
    sample: str,
    year: str,
    era: str,
    primary_dataset: str,
    is_mc: bool = False,
    nevents: str = "0",
    nano_version: int = 15,
    extra_metadata: dict[str, str | int | bool] | None = None,
) -> dict:
    """Build a PocketCoffea dataset-definition dictionary."""
    metadata: dict[str, str | int | bool] = {
        "sample": sample,
        "year": year,
        "era": era,
        "isMC": str(is_mc),
        "primaryDataset": primary_dataset,
        "nevents": str(nevents),
        "nano_version": nano_version,
    }
    if extra_metadata:
        metadata.update(extra_metadata)

    return {
        dataset_name: {
            "files": files,
            "metadata": metadata,
        }
    }


def write_dataset_definition(dataset: dict, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(dataset, indent=2) + "\n")
