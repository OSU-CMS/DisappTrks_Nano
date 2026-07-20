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
    nano_version: int


ERA_GROUPS = (
    EraGroup("2022CD", (("2022", "C"), ("2022", "D")), "2022_preEE", "CD", 12),
    EraGroup(
        "2022EFG",
        (("2022", "E"), ("2022", "F"), ("2022", "G")),
        "2022_postEE",
        "EFG",
        12,
    ),
    EraGroup("2023C", (("2023", "C"),), "2023_preBPix", "C", 12),
    EraGroup("2023D", (("2023", "D"),), "2023_postBPix", "D", 12),
    EraGroup(
        "2024",
        tuple(("2024", era) for era in "ABCDEFGHI"),
        "2024",
        "all",
        15,
    ),
    EraGroup(
        "2025",
        tuple(("2025", era) for era in "ABCDEFGHI"),
        "2025",
        "all",
        15,
    ),
    EraGroup(
        "2026",
        tuple(("2026", era) for era in ("ABCD")),
        "2026",
        "all",
        15,
    )
)

ERA_GROUP_BY_LABEL = {group.label: group for group in ERA_GROUPS}
YEAR_ERA_TO_GROUP = {
    year_era: group.label for group in ERA_GROUPS for year_era in group.years_eras
}
PRIMARY_DATASETS = ("Muon", "EGamma", "JetMET")
ALLOWED_DEV_DIRS = (
    "EGamma0",
    "EGamma1",
    "EGamma2",
    "EGamma22",
    "EGamma3",
    "EGamma4",
    "EGamma5",
    "JetMET",
    "JetMET0",
    "JetMET1",
    "JetMET2",
    "JetMET3",
    "Muon",
    "Muon0",
    "Muon1",
    "Muon3",
)
ALLOWED_PROD_DIRS = (
    "JetMET_Run2022C",
    "JetMET_Run2022D",
    "JetMET_Run2022E",
    "JetMET_Run2022F",
    "JetMET_Run2022G",
    "JetMET_Run2023C",
    "JetMET_Run2023C_v1",
    "JetMET_Run2023C_v2",
    "JetMET_Run2023C_v3",
    "JetMET_Run2023C_v4",
    "JetMET_Run2023D",
    "JetMET_Run2023D_v1",
    "JetMET_Run2023D_v2",
    "JetMET_Run2024C",
    "JetMET_Run2024D",
    "JetMET_Run2024E",
    "JetMET_Run2024F",
    "JetMET_Run2024G",
    "JetMET_Run2024H",
    "JetMET_Run2024I",
    "JetMET0_Run2023C_v1",
    "JetMET0_Run2023C_v2",
    "JetMET0_Run2023C_v3",
    "JetMET0_Run2023C_v4",
    "JetMET0_Run2023D_v1",
    "JetMET0_Run2023D_v2",
    "JetMET0_Run2024C",
    "JetMET0_Run2024D",
    "JetMET0_Run2024E",
    "JetMET0_Run2024F",
    "JetMET0_Run2024G",
    "JetMET0_Run2024H",
    "JetMET0_Run2024I",
    "JetMET1_Run2023C_v1",
    "JetMET1_Run2023C_v2",
    "JetMET1_Run2023C_v3",
    "JetMET1_Run2023D_v1",
    "JetMET1_Run2023D_v2",
    "JetMET1_Run2024C",
    "JetMET1_Run2024D",
    "JetMET1_Run2024E",
    "JetMET1_Run2024F",
    "Muon0_Run2023C_v1",
    "Muon0_Run2023C_v2",
    "Muon0_Run2023C_v3",
    "Muon0_Run2023C_v4",
    "Muon0_Run2023D_v1",
    "Muon0_Run2023D_v2",
    "Muon0_Run2024C",
    "Muon0_Run2024D",
    "Muon0_Run2024E",
    "Muon0_Run2024F",
    "Muon0_Run2024G",
    "Muon0_Run2024H",
    "Muon0_Run2024I",
    "Muon0_Run2024I_v2",
    "Muon1_Run2023C_v1",
    "Muon1_Run2023C_v2",
    "Muon1_Run2023C_v3",
    "Muon1_Run2023D_v1",
    "Muon1_Run2023D_v2",
    "Muon1_Run2024C",
    "Muon1_Run2024D",
    "Muon1_Run2024E",
    "Muon1_Run2024F",
    "Muon_Run2022C",
    "Muon_Run2022D",
    "Muon_Run2022F",
    "Muon_Run2022G",
)
RUN_RE = re.compile(r"Run(20\d{2})([A-Z])")
PROD_DATASET_DIR_RE = re.compile(
    r"^(?P<primary>Muon\d*|EGamma\d*|JetMET\d*)_Run(?P<year>20\d{2})(?P<era>[A-Z])(?:_v(?P<version>\d+))?$"
)


def primary_dataset_from_path(path: str) -> str | None:
    """Infer Muon/EGamma/JetMET from an OSUNano EOS path."""
    for part in Path(path).parts:
        if part.startswith("Muon"):
            return "Muon"
        if part.startswith("EGamma"):
            return "EGamma"
        if part.startswith("JetMET"):
            return "JetMET"
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


def osunano_area_and_top_dir(path: str) -> tuple[str, str] | None:
    """Return ``(dev|prod, top_dir)`` for files under the OSUNano EOS areas."""
    parts = Path(path).parts
    for area in ("dev", "prod"):
        if area in parts:
            index = parts.index(area)
            if index + 1 < len(parts):
                return area, parts[index + 1]
    return None


def is_allowed_osunano_path(
    path: str,
    *,
    allowed_dev_dirs: tuple[str, ...] = ALLOWED_DEV_DIRS,
    allowed_prod_dirs: tuple[str, ...] = ALLOWED_PROD_DIRS,
) -> bool:
    """Keep only files from the explicitly listed case-sensitive EOS dirs."""
    area_top = osunano_area_and_top_dir(path)
    if area_top is None:
        return False
    area, top_dir = area_top
    if top_dir.startswith("JetMET"):
        return True
    if area == "dev":
        return top_dir in allowed_dev_dirs
    if area == "prod":
        return top_dir in allowed_prod_dirs
    return False


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
    prod_version_policy: str = "all",
) -> dict[tuple[str, str], list[str]]:
    """Group OSUNano ROOT files by ``(primary_dataset, era_group_label)``."""
    if prod_version_policy not in ("latest", "all"):
        raise ValueError("prod_version_policy must be 'latest' or 'all'")
    files = [path for path in files if is_allowed_osunano_path(path)]
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
    nano_version: int | None = None,
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
                nano_version=group.nano_version if nano_version is None else nano_version,
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
