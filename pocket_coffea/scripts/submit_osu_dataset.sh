#!/bin/bash
set -euo pipefail

if [ "$#" -lt 4 ]; then
    cat <<'USAGE' >&2
Usage:
  scripts/submit_osu_dataset.sh DATASET_JSON TAG SAMPLE YEAR

Example:
  FILES_PER_JOB=5 CHUNKSIZE=50000 CATEGORY_MODE=muon_pveto \
    scripts/submit_osu_dataset.sh datasets/osu_2023C_muon.json 2023C_muon_pveto DATA_Muon 2023_preBPix

Environment overrides:
  FILES_PER_JOB    default: 5
  CHUNKSIZE        default: 50000
  JOB_TIMEOUT      default: 2h
  MAX_JOBS         default: all jobs
  REQUEST_MEMORY   default: 4GB
  REQUEST_DISK     default: 4GB
  CATEGORY_MODE    default: muon_pveto
  TRANSFER_ROOT_FILES default: 0; set to 1 when worker nodes cannot see dataset paths
  X509_USER_PROXY  optional; only transferred if it points to an existing file
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
TRANSFER_ROOT_FILES="${TRANSFER_ROOT_FILES:-0}"
X509_PROXY="${X509_USER_PROXY:-}"

if [ "${TRANSFER_ROOT_FILES}" = "1" ] && [ "${FILES_PER_JOB}" -ne 1 ]; then
    cat <<'MSG' >&2
TRANSFER_ROOT_FILES=1 currently requires FILES_PER_JOB=1.

This avoids Condor sandbox filename collisions from inputs like many different
directories containing nano_1.root.  Start with FILES_PER_JOB=1; we can add
staging-directory support later if needed.
MSG
    exit 1
fi

if [ ! -f "${DATASET_JSON}" ]; then
    echo "Dataset JSON not found: ${DATASET_JSON}" >&2
    exit 1
fi

if [ ! -d "python_env" ]; then
    cat <<'MSG' >&2
Missing pocket_coffea/python_env.

Build a relocatable Python target directory before submitting, e.g. from the
repository root:

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
    echo "Warning: data/golden_jsons not found; data jobs may fail if /cvmfs golden JSONs are unavailable." >&2
fi
if [ -d "data/jet_veto_maps" ]; then
    TRANSFER_INPUTS="${TRANSFER_INPUTS},data/jet_veto_maps"
else
    echo "Warning: data/jet_veto_maps not found; data jobs may fail if /cvmfs JME jet-veto-map payloads are unavailable." >&2
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
echo "  transfer ROOT: ${TRANSFER_ROOT_FILES}"
if [ "${X509_BASENAME}" != "__none__" ]; then
    echo "  proxy:         ${X509_PROXY_ABS}"
else
    echo "  proxy:         none"
fi

if [ "${TRANSFER_ROOT_FILES}" = "1" ]; then
    QUEUE_FILE="logs/${TAG}/root_transfer_queue.txt"
    SUBMIT_FILE="logs/${TAG}/submit_root_transfer.jdl"
    python3 - "${DATASET_JSON_ABS}" "${FILES_PER_JOB}" "${NJOBS}" > "${QUEUE_FILE}" <<'PY'
import json
import sys

dataset_json, files_per_job, n_jobs = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
with open(dataset_json) as handle:
    dataset = json.load(handle)

files = []
for definition in dataset.values():
    files.extend(definition.get("files", []))

for job_id in range(n_jobs):
    start = job_id * files_per_job
    selected = files[start:start + files_per_job]
    if not selected:
        continue
    print(job_id, ",".join(selected))
PY

    cat > "${SUBMIT_FILE}" <<EOF
universe = vanilla
executable = scripts/run_osu_pocket_coffea_job.sh
arguments = \$(job_id) ${DATASET_BASENAME} ${FILES_PER_JOB} ${CHUNKSIZE} ${SAMPLE} ${YEAR} ${TAG} ${X509_BASENAME} ${JOB_TIMEOUT} ${CATEGORY_MODE} 1

output = logs/${TAG}/job_\$(ClusterId).\$(job_id).out
error  = logs/${TAG}/job_\$(ClusterId).\$(job_id).err
log    = logs/${TAG}/job_\$(ClusterId).log

request_cpus = 1
request_memory = ${REQUEST_MEMORY}
request_disk = ${REQUEST_DISK}

should_transfer_files = YES
when_to_transfer_output = ON_EXIT
transfer_input_files = ${TRANSFER_INPUTS},\$(root_inputs)
transfer_output_files = analysis_output

on_exit_remove = (ExitBySignal == False) && (ExitCode == 0)
max_retries = 3
requirements = Machine =!= LastRemoteHost

queue job_id,root_inputs from ${QUEUE_FILE}
EOF

    condor_submit "${SUBMIT_FILE}"
else
    condor_submit \
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
        -append "n_jobs=${NJOBS}" \
        scripts/submit_osu_pocket_coffea.jdl
fi
