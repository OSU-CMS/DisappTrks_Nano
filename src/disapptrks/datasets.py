"""Helpers for building PocketCoffea dataset JSON files."""

from __future__ import annotations

import json
import re
import subprocess
from dataclasses import dataclass
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


@dataclass(frozen=True)
class EraGroup:
    label: str
    years_eras: tuple[tuple[str, str], ...]
    metadata_year: str
    metadata_era: str


ERA_GROUPS = (
    EraGroup("2022CD", (("2022", "C"), ("2022", "D")), "2022_preEE", "CD"),
    EraGroup(
        "2022EFG",
        (("2022", "E"), ("2022", "F"), ("2022", "G")),
        "2022_postEE",
        "EFG",
    ),
    EraGroup("2023C", (("2023", "C"),), "2023_preBPix", "C"),
    EraGroup("2023D", (("2023", "D"),), "2023_postBPix", "D"),
    EraGroup(
        "2024",
        tuple(("2024", era) for era in "ABCDEFGHI"),
        "2024",
        "all",
    ),
    EraGroup(
        "2025",
        tuple(("2025", era) for era in "ABCDEFGHI"),
        "2025",
        "all",
    ),
)

ERA_GROUP_BY_LABEL = {group.label: group for group in ERA_GROUPS}
YEAR_ERA_TO_GROUP = {
    year_era: group.label for group in ERA_GROUPS for year_era in group.years_eras
}
PRIMARY_DATASETS = ("Muon", "EGamma")
RUN_RE = re.compile(r"Run(20\d{2})([A-Z])")
PROD_DATASET_DIR_RE = re.compile(
    r"^(?P<primary>Muon\d*|EGamma\d*)_Run(?P<year>20\d{2})(?P<era>[A-Z])(?:_v(?P<version>\d+))?$"
)


def primary_dataset_from_path(path: str) -> str | None:
    """Infer Muon/EGamma from an OSUNano EOS path."""
    for part in Path(path).parts:
        if part.startswith("Muon"):
            return "Muon"
        if part.startswith("EGamma"):
            return "EGamma"
    return None


def run_year_era_from_path(path: str) -> tuple[str, str] | None:
    """Return the first ``(year, era)`` found in a RunYYYYX path or filename."""
    match = RUN_RE.search(path)
    if not match:
        return None
    return match.group(1), match.group(2)


def era_group_label_from_path(path: str) -> str | None:
    year_era = run_year_era_from_path(path)
    if year_era is None:
        return None
    return YEAR_ERA_TO_GROUP.get(year_era)


def _prod_dataset_component(path: str) -> str | None:
    parts = Path(path).parts
    if "prod" not in parts:
        return None
    for part in parts:
        if PROD_DATASET_DIR_RE.match(part):
            return part
    return None


def _prod_dataset_base_and_version(component: str) -> tuple[str, int]:
    match = PROD_DATASET_DIR_RE.match(component)
    if match is None:
        return component, 1
    version = int(match.group("version") or 1)
    base = re.sub(r"_v\d+$", "", component)
    return base, version


def filter_latest_prod_versions(files: list[str]) -> list[str]:
    """Drop older ``*_vN`` prod directories while keeping dev files untouched."""
    latest: dict[str, tuple[int, str]] = {}
    for path in files:
        component = _prod_dataset_component(path)
        if component is None:
            continue
        base, version = _prod_dataset_base_and_version(component)
        if base not in latest or version > latest[base][0]:
            latest[base] = (version, component)

    if not latest:
        return files

    allowed_components = {component for _, component in latest.values()}
    filtered = []
    for path in files:
        component = _prod_dataset_component(path)
        if component is None or component in allowed_components:
            filtered.append(path)
    return filtered


def group_osunano_files(
    files: list[str],
    *,
    primary_datasets: tuple[str, ...] = PRIMARY_DATASETS,
    group_labels: tuple[str, ...] = tuple(group.label for group in ERA_GROUPS),
    prod_version_policy: str = "latest",
) -> dict[tuple[str, str], list[str]]:
    """Group OSUNano ROOT files by ``(primary_dataset, era_group_label)``."""
    if prod_version_policy not in ("latest", "all"):
        raise ValueError("prod_version_policy must be 'latest' or 'all'")
    if prod_version_policy == "latest":
        files = filter_latest_prod_versions(files)

    grouped = {
        (primary, label): []
        for primary in primary_datasets
        for label in group_labels
    }
    for path in files:
        primary = primary_dataset_from_path(path)
        label = era_group_label_from_path(path)
        if primary in primary_datasets and label in group_labels:
            grouped[(primary, label)].append(path)

    return {key: sorted(value) for key, value in grouped.items()}


def write_filelist(files: list[str], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("".join(f"{path}\n" for path in sorted(files)))


def scan_eos_bases_for_root_files(
    eos_paths: list[str],
    *,
    xrootd: str = "root://cmseos.fnal.gov",
) -> list[str]:
    files: list[str] = []
    for eos_path in eos_paths:
        files.extend(list_eos_root_files(eos_path, xrootd=xrootd, recursive=True))
    return sorted(set(files))


def write_grouped_filelists(
    grouped: dict[tuple[str, str], list[str]],
    *,
    output_dir: Path,
    dataset_json_dir: Path | None = None,
    nano_version: int = 15,
) -> dict[str, dict[str, Path | int | str]]:
    """Write grouped filelists and optionally matching PocketCoffea JSONs."""
    outputs: dict[str, dict[str, Path | int | str]] = {}
    for (primary, label), files in sorted(grouped.items()):
        key = f"{primary}_{label}"
        filelist_path = output_dir / f"{primary}_{label}.txt"
        write_filelist(files, filelist_path)

        group = ERA_GROUP_BY_LABEL[label]
        entry: dict[str, Path | int | str] = {
            "filelist": filelist_path,
            "n_files": len(files),
        }

        if dataset_json_dir is not None:
            dataset_name = f"Run{label}_{primary}_OSUNano_EOS"
            sample = f"DATA_{primary}"
            dataset_path = dataset_json_dir / f"eos_{label}_{primary}.json"
            dataset = build_dataset_definition(
                dataset_name=dataset_name,
                files=files,
                sample=sample,
                year=group.metadata_year,
                era=group.metadata_era,
                primary_dataset=primary,
                is_mc=False,
                nevents="0",
                nano_version=nano_version,
            )
            write_dataset_definition(dataset, dataset_path)
            entry["dataset_json"] = dataset_path
            entry["dataset_name"] = dataset_name
            entry["sample"] = sample
            entry["year"] = group.metadata_year
            entry["era"] = group.metadata_era

        outputs[key] = entry
    return outputs


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
