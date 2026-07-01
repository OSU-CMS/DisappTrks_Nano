#!/bin/bash
set -euo pipefail

JOBID="$1"
DATASET_JSON="$2"
FILES_PER_JOB="$3"
CHUNKSIZE="$4"
SAMPLE="$5"
YEAR="$6"
TAG="$7"
X509_BASENAME="${8:-}"

export XRD_RUNFORKHANDLER=1
export MALLOC_TRIM_THRESHOLD_=0

if [ -n "${X509_BASENAME}" ] && [ -f "${X509_BASENAME}" ]; then
    export X509_USER_PROXY="$PWD/${X509_BASENAME}"
elif [ -f "x509up" ]; then
    export X509_USER_PROXY="$PWD/x509up"
elif [ -n "${X509_USER_PROXY:-}" ] && [ -f "$X509_USER_PROXY" ]; then
    export X509_USER_PROXY="$X509_USER_PROXY"
fi

echo "Starting LPC PocketCoffea job ${JOBID} on $(hostname)"
echo "Sandbox: $PWD"
echo "Dataset JSON: ${DATASET_JSON}"
echo "Files per job: ${FILES_PER_JOB}"
echo "Chunksize: ${CHUNKSIZE}"
echo "Sample/year/tag: ${SAMPLE} ${YEAR} ${TAG}"
echo "X509_USER_PROXY: ${X509_USER_PROXY:-unset}"
python3 --version

# Create the transfer-output directory before doing anything fragile.  If the
# job fails during environment setup or imports, Condor can still transfer this
# directory and leave the real failure in the .err/.out logs instead of holding
# the job because ``analysis_output`` does not exist.
mkdir -p "analysis_output/${TAG}"

export PYTHONPATH="$PWD/python_env:$PWD/src:$PWD:${PYTHONPATH:-}"
export PATH="$PWD/python_env/bin:${PATH}"

python3 - <<'PY'
import importlib
for name in ("awkward", "uproot", "coffea", "pocket_coffea", "disapptrks"):
    module = importlib.import_module(name)
    version = getattr(module, "__version__", "unknown")
    print(f"import {name}: {version}")
PY

JOB_DATASET_JSON="job_dataset_${JOBID}.json"
python3 make_lpc_job_dataset.py \
    --input "${DATASET_JSON}" \
    --output "${JOB_DATASET_JSON}" \
    --job-id "${JOBID}" \
    --files-per-job "${FILES_PER_JOB}"

export DISAPPTRKS_DATASET_JSON="$PWD/${JOB_DATASET_JSON}"
export DISAPPTRKS_DATASET_SAMPLE="${SAMPLE}"
export DISAPPTRKS_DATASET_YEAR="${YEAR}"
unset DISAPPTRKS_ENABLE_SEARCH_DIAGNOSTICS

JOB_OUTPUT="job_output_${JOBID}"
python3 -m pocket_coffea.scripts.runner \
    --cfg config.py \
    --process-separately \
    --executor iterative \
    --chunksize "${CHUNKSIZE}" \
    -o "${JOB_OUTPUT}"

shopt -s nullglob
outputs=("${JOB_OUTPUT}"/*.coffea)
if [ "${#outputs[@]}" -eq 0 ]; then
    echo "No .coffea outputs found in ${JOB_OUTPUT}" >&2
    exit 1
fi

for output in "${outputs[@]}"; do
    base="$(basename "${output}" .coffea)"
    cp "${output}" "analysis_output/${TAG}/${base}_job_${JOBID}.coffea"
    echo "Saved analysis_output/${TAG}/${base}_job_${JOBID}.coffea"
done

echo "Done"
