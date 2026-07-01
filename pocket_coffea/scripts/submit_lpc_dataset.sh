#!/bin/bash
set -euo pipefail

if [ "$#" -lt 4 ]; then
    cat <<'USAGE' >&2
Usage:
  scripts/submit_lpc_dataset.sh DATASET_JSON TAG SAMPLE YEAR

Example:
  FILES_PER_JOB=5 CHUNKSIZE=50000 \
    scripts/submit_lpc_dataset.sh datasets/eos_2023C_muon.json 2023C_muon_pveto DATA_Muon 2023

Environment overrides:
  FILES_PER_JOB    default: 5
  CHUNKSIZE        default: 50000
  REQUEST_MEMORY   default: 4GB
  REQUEST_DISK     default: 2GB
USAGE
    exit 2
fi

DATASET_JSON="$1"
TAG="$2"
SAMPLE="$3"
YEAR="$4"

FILES_PER_JOB="${FILES_PER_JOB:-5}"
CHUNKSIZE="${CHUNKSIZE:-50000}"
REQUEST_MEMORY="${REQUEST_MEMORY:-4GB}"
REQUEST_DISK="${REQUEST_DISK:-2GB}"
X509_PROXY="${X509_USER_PROXY:-}"

if [ -z "${X509_PROXY}" ] && command -v voms-proxy-info >/dev/null 2>&1; then
    X509_PROXY="$(voms-proxy-info -path 2>/dev/null || true)"
fi

if [ -z "${X509_PROXY}" ]; then
    X509_PROXY="${HOME}/x509up_u$(id -u)"
fi

if [ ! -f "${DATASET_JSON}" ]; then
    echo "Dataset JSON not found: ${DATASET_JSON}" >&2
    exit 1
fi

if [ ! -d "python_env" ]; then
    cat <<'MSG' >&2
Missing pocket_coffea/python_env.

Build or copy a relocatable Python target directory before submitting, e.g.
from the repository root:

  python3 -m pip install --target pocket_coffea/python_env '.[analysis,pocket-coffea]'

MSG
    exit 1
fi

if [ ! -d "../src/disapptrks" ]; then
    echo "Run this from DisappTrks_Nano/pocket_coffea so ../src/disapptrks exists." >&2
    exit 1
fi

if [ ! -f "${X509_PROXY}" ]; then
    cat <<MSG >&2
Could not find an X509 proxy to transfer:
  ${X509_PROXY}

Create one with voms-proxy-init, or set X509_USER_PROXY to the proxy path.
MSG
    exit 1
fi

NFILES="$(python3 - "${DATASET_JSON}" <<'PY'
import json
import sys
with open(sys.argv[1]) as handle:
    dataset = json.load(handle)
print(sum(len(definition.get("files", [])) for definition in dataset.values()))
PY
)"

if [ "${NFILES}" -le 0 ]; then
    echo "No files found in ${DATASET_JSON}" >&2
    exit 1
fi

NJOBS=$(( (NFILES + FILES_PER_JOB - 1) / FILES_PER_JOB ))
DATASET_JSON_ABS="$(python3 - "${DATASET_JSON}" <<'PY'
from pathlib import Path
import sys
print(Path(sys.argv[1]).resolve())
PY
)"
DATASET_BASENAME="$(basename "${DATASET_JSON_ABS}")"
X509_PROXY_ABS="$(python3 - "${X509_PROXY}" <<'PY'
from pathlib import Path
import sys
print(Path(sys.argv[1]).resolve())
PY
)"
X509_BASENAME="$(basename "${X509_PROXY_ABS}")"

mkdir -p "logs/${TAG}" "analysis_output/${TAG}"

echo "Submitting ${NJOBS} job(s) for ${NFILES} file(s)"
echo "  dataset:       ${DATASET_JSON_ABS}"
echo "  tag:           ${TAG}"
echo "  sample/year:   ${SAMPLE} / ${YEAR}"
echo "  files/job:     ${FILES_PER_JOB}"
echo "  chunksize:     ${CHUNKSIZE}"
echo "  memory/disk:   ${REQUEST_MEMORY} / ${REQUEST_DISK}"
echo "  proxy:         ${X509_PROXY_ABS}"

condor_submit \
    -append "dataset_json=${DATASET_JSON_ABS}" \
    -append "dataset_basename=${DATASET_BASENAME}" \
    -append "x509_proxy=${X509_PROXY_ABS}" \
    -append "x509_basename=${X509_BASENAME}" \
    -append "tag=${TAG}" \
    -append "sample=${SAMPLE}" \
    -append "year=${YEAR}" \
    -append "files_per_job=${FILES_PER_JOB}" \
    -append "chunksize=${CHUNKSIZE}" \
    -append "request_memory=${REQUEST_MEMORY}" \
    -append "request_disk=${REQUEST_DISK}" \
    -append "n_jobs=${NJOBS}" \
    scripts/submit_lpc_pocket_coffea.jdl
