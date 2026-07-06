#!/bin/bash
set -euo pipefail

if [ "$#" -lt 4 ]; then
    cat <<'USAGE' >&2
Usage:
  scripts/submit_lxplus_dataset.sh DATASET_JSON TAG SAMPLE YEAR

Example:
  FILES_PER_JOB=5 CHUNKSIZE=50000 CATEGORY_MODE=muon_pveto \
    scripts/submit_lxplus_dataset.sh datasets/eos_2023C_muon.json 2023C_muon_pveto DATA_Muon 2023_preBPix

Environment overrides:
  FILES_PER_JOB    default: 5
  CHUNKSIZE        default: 50000
  JOB_TIMEOUT      default: 2h
  MAX_JOBS         default: all jobs
  REQUEST_MEMORY   default: 4GB
  REQUEST_DISK     default: 4GB
  CATEGORY_MODE    default: muon_pveto
  JOB_FLAVOUR      default: workday
  WANT_OS          default: el9
  OUTPUT_DEST      optional; copy .coffea files there as well as Condor transfer
  FILE_REWRITE_FROM default: root://cmseosmgm01.fnal.gov:1094//
  FILE_REWRITE_TO   default: root://cmseos.fnal.gov//
  DRY_RUN          default: 0; print resolved settings without condor_submit
  X509_USER_PROXY  optional; transferred if it points to an existing file
USAGE
    exit 2
fi

DATASET_JSON="$1"
TAG="$2"
SAMPLE="$3"
YEAR="$4"

FILES_PER_JOB="${FILES_PER_JOB:-5}"
CHUNKSIZE="${CHUNKSIZE:-50000}"
JOB_TIMEOUT="${JOB_TIMEOUT:-2h}"
MAX_JOBS="${MAX_JOBS:-}"
REQUEST_MEMORY="${REQUEST_MEMORY:-4GB}"
REQUEST_DISK="${REQUEST_DISK:-4GB}"
CATEGORY_MODE="${CATEGORY_MODE:-muon_pveto}"
JOB_FLAVOUR="${JOB_FLAVOUR:-workday}"
WANT_OS="${WANT_OS:-el9}"
OUTPUT_DEST="${OUTPUT_DEST:-}"
FILE_REWRITE_FROM="${FILE_REWRITE_FROM:-root://cmseosmgm01.fnal.gov:1094//}"
FILE_REWRITE_TO="${FILE_REWRITE_TO:-root://cmseos.fnal.gov//}"
DRY_RUN="${DRY_RUN:-0}"
X509_PROXY="${X509_USER_PROXY:-}"

submit_with_retries() {
    local attempt
    local retries="${SUBMIT_RETRIES:-5}"

    for attempt in $(seq 1 "${retries}"); do
        echo "condor_submit scripts/submit_lxplus_pocket_coffea.jdl (attempt ${attempt}/${retries})"
        if condor_submit \
            -append "dataset_basename=${DATASET_BASENAME}" \
            -append "transfer_inputs=${TRANSFER_INPUTS}" \
            -append "x509_basename=${X509_BASENAME}" \
            -append "tag=${TAG}" \
            -append "sample=${SAMPLE}" \
            -append "year=${YEAR}" \
            -append "files_per_job=${FILES_PER_JOB}" \
            -append "chunksize=${CHUNKSIZE}" \
            -append "job_timeout=${JOB_TIMEOUT}" \
            -append "request_memory=${REQUEST_MEMORY}" \
            -append "request_disk=${REQUEST_DISK}" \
            -append "category_mode=${CATEGORY_MODE}" \
            -append "job_flavour=${JOB_FLAVOUR}" \
            -append "want_os=${WANT_OS}" \
            -append "output_dest=${OUTPUT_DEST:-__none__}" \
            -append "file_rewrite_from=${FILE_REWRITE_FROM:-__none__}" \
            -append "file_rewrite_to=${FILE_REWRITE_TO:-__none__}" \
            -append "n_jobs=${NJOBS}" \
            scripts/submit_lxplus_pocket_coffea.jdl
        then
            return 0
        fi
        if [ "${attempt}" -lt "${retries}" ]; then
            sleep 10
        fi
    done

    echo "ERROR: failed to submit scripts/submit_lxplus_pocket_coffea.jdl after ${retries} attempts" >&2
    return 1
}

if [ ! -f "${DATASET_JSON}" ]; then
    echo "Dataset JSON not found: ${DATASET_JSON}" >&2
    exit 1
fi

if [ ! -d "python_env" ]; then
    cat <<'MSG' >&2
Missing pocket_coffea/python_env.

Build a relocatable Python target directory on lxplus before submitting, e.g.
from the repository root:

  python3 -m pip install --target pocket_coffea/python_env '.[analysis,pocket-coffea]'

MSG
    exit 1
fi

if [ ! -d "../src/disapptrks" ]; then
    echo "Run this from DisappTrks_Nano/pocket_coffea so ../src/disapptrks exists." >&2
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
if [ -n "${MAX_JOBS}" ]; then
    if [ "${MAX_JOBS}" -le 0 ]; then
        echo "MAX_JOBS must be positive if set: ${MAX_JOBS}" >&2
        exit 1
    fi
    if [ "${MAX_JOBS}" -lt "${NJOBS}" ]; then
        NJOBS="${MAX_JOBS}"
    fi
fi

DATASET_JSON_ABS="$(python3 - "${DATASET_JSON}" <<'PY'
from pathlib import Path
import sys
print(Path(sys.argv[1]).resolve())
PY
)"
DATASET_BASENAME="$(basename "${DATASET_JSON_ABS}")"

TRANSFER_INPUTS="config.py,cuts.py,workflow.py,../src,python_env,${DATASET_JSON_ABS},scripts/make_lpc_job_dataset.py"
if [ -d "data/golden_jsons" ]; then
    TRANSFER_INPUTS="${TRANSFER_INPUTS},data/golden_jsons"
else
    echo "Warning: data/golden_jsons not found; data jobs may rely on CVMFS golden JSONs." >&2
fi
if [ -d "data/jet_veto_maps" ]; then
    TRANSFER_INPUTS="${TRANSFER_INPUTS},data/jet_veto_maps"
else
    echo "Warning: data/jet_veto_maps not found; data jobs may rely on CVMFS JME jet-veto-map payloads." >&2
fi

X509_BASENAME="__none__"
if [ -n "${X509_PROXY}" ] && [ -f "${X509_PROXY}" ]; then
    X509_PROXY_ABS="$(python3 - "${X509_PROXY}" <<'PY'
from pathlib import Path
import sys
print(Path(sys.argv[1]).resolve())
PY
)"
    X509_BASENAME="$(basename "${X509_PROXY_ABS}")"
    TRANSFER_INPUTS="${TRANSFER_INPUTS},${X509_PROXY_ABS}"
fi

mkdir -p "logs/${TAG}" "analysis_output/${TAG}"

echo "Submitting ${NJOBS} job(s) for ${NFILES} file(s)"
echo "  dataset:       ${DATASET_JSON_ABS}"
echo "  tag:           ${TAG}"
echo "  sample/year:   ${SAMPLE} / ${YEAR}"
echo "  files/job:     ${FILES_PER_JOB}"
if [ -n "${MAX_JOBS}" ]; then
    echo "  max jobs:      ${MAX_JOBS}"
fi
echo "  chunksize:     ${CHUNKSIZE}"
echo "  job timeout:   ${JOB_TIMEOUT}"
echo "  memory/disk:   ${REQUEST_MEMORY} / ${REQUEST_DISK}"
echo "  category mode: ${CATEGORY_MODE}"
echo "  job flavour:   ${JOB_FLAVOUR}"
echo "  want OS:       ${WANT_OS}"
echo "  output dest:   ${OUTPUT_DEST:-Condor transfer only}"
if [ -n "${FILE_REWRITE_FROM}" ]; then
    echo "  file rewrite:  ${FILE_REWRITE_FROM} -> ${FILE_REWRITE_TO}"
else
    echo "  file rewrite:  none"
fi
echo "  dry run:       ${DRY_RUN}"
if [ "${X509_BASENAME}" != "__none__" ]; then
    echo "  proxy:         ${X509_PROXY_ABS}"
else
    echo "  proxy:         none"
fi

if [ "${DRY_RUN}" = "1" ]; then
    echo
    echo "DRY_RUN=1: not submitting. Resolved transfer inputs:"
    echo "  ${TRANSFER_INPUTS}"
    exit 0
fi

submit_with_retries
