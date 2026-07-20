#!/usr/bin/env python3
"""Download golden JSON files used by the Run-3 data preselections."""

from __future__ import annotations

import argparse
import urllib.request
from pathlib import Path


GOLDEN_JSON_URLS = {
    "Cert_Collisions2022_355100_362760_Golden.json": (
        "https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions22/"
        "Cert_Collisions2022_355100_362760_Golden.json"
    ),
    "Cert_Collisions2023_366442_370790_Golden.json": (
        "https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions23/"
        "Cert_Collisions2023_366442_370790_Golden.json"
    ),
    "Cert_Collisions2024_378981_386951_Golden.json": (
        "https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions24/"
        "Cert_Collisions2024_378981_386951_Golden.json"
    ),
    "Cert_Collisions2025_391658_398903_Golden.json": (
        "https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions25/"
        "Cert_Collisions2025_391658_398903_Golden.json"
    ),
    "Collisions26_MLEnhancedGolden_Latest.json": (
        "https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions26/"
        "Collisions26_MLEnhancedGolden_Latest.json"
    ),
}


def main() -> int:
    parser = argparse.ArgumentParser(description="Download CMS golden JSON files.")
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        default=Path("data/golden_jsons"),
        help="Destination directory. Default: data/golden_jsons.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Redownload files that already exist.",
    )
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    for filename, url in GOLDEN_JSON_URLS.items():
        output = args.output_dir / filename
        if output.exists() and not args.overwrite:
            print(f"skip existing {output}")
            continue
        print(f"download {url} -> {output}")
        urllib.request.urlretrieve(url, output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
