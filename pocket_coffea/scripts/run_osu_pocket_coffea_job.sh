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
JOB_TIMEOUT="${9:-2h}"
CATEGORY_MODE="${10:-muon_pveto}"

export XRD_RUNFORKHANDLER=1
export MALLOC_TRIM_THRESHOLD_=0

if [ -n "${X509_BASENAME}" ] && [ "${X509_BASENAME}" != "__none__" ] && [ -f "${X509_BASENAME}" ]; then
    export X509_USER_PROXY="$PWD/${X509_BASENAME}"
elif [ -n "${X509_USER_PROXY:-}" ] && [ -f "${X509_USER_PROXY}" ]; then
    export X509_USER_PROXY="$X509_USER_PROXY"
fi

mkdir -p "analysis_output/${TAG}"
exec > >(tee -a "analysis_output/${TAG}/job_${JOBID}.debug.log") 2>&1

on_error() {
    local status="$?"
    echo "Job ${JOBID} failed with exit code ${status}"
    echo "PWD: $PWD"
    echo "Directory listing:"
    ls -la
    echo "analysis_output listing:"
    find analysis_output -maxdepth 3 -type f -print || true
    exit "${status}"
}
trap on_error ERR

echo "Starting OSU PocketCoffea job ${JOBID} on $(hostname)"
echo "Sandbox: $PWD"
echo "Dataset JSON: ${DATASET_JSON}"
echo "Files per job: ${FILES_PER_JOB}"
echo "Chunksize: ${CHUNKSIZE}"
echo "Job timeout: ${JOB_TIMEOUT}"
echo "Sample/year/tag: ${SAMPLE} ${YEAR} ${TAG}"
echo "Category mode: ${CATEGORY_MODE}"
echo "X509_USER_PROXY: ${X509_USER_PROXY:-unset}"
python3 --version
echo "Initial sandbox contents:"
ls -la

export PYTHONPATH="$PWD/src:$PWD:$PWD/python_env:${PYTHONPATH:-}"
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
    --files-per-job "${FILES_PER_JOB}" \
    --fallback-events-per-file "${CHUNKSIZE}"

export DISAPPTRKS_DATASET_JSON="$PWD/${JOB_DATASET_JSON}"
export DISAPPTRKS_DATASET_SAMPLE="${SAMPLE}"
export DISAPPTRKS_DATASET_YEAR="${YEAR}"
export DISAPPTRKS_CATEGORY_MODE="${CATEGORY_MODE}"
if [ -d "$PWD/data/golden_jsons" ]; then
    export DISAPPTRKS_GOLDEN_JSON_DIR="$PWD/data/golden_jsons"
fi
unset DISAPPTRKS_ENABLE_SEARCH_DIAGNOSTICS

cat > osu_inner_run_options.yaml <<'YAML'
ignore-grid-certificate: true
YAML

JOB_OUTPUT="job_output_${JOBID}"
runner_cmd=(
    python3 -m pocket_coffea.scripts.runner
    --cfg config.py
    --process-separately
    --executor iterative
    --chunksize "${CHUNKSIZE}"
    --custom-run-options osu_inner_run_options.yaml
    -o "${JOB_OUTPUT}"
)

if command -v timeout >/dev/null 2>&1; then
    timeout --preserve-status "${JOB_TIMEOUT}" "${runner_cmd[@]}"
else
    echo "WARNING: timeout command not found; running without a walltime guard"
    "${runner_cmd[@]}"
fi

echo "PocketCoffea output tree:"
find "${JOB_OUTPUT}" -maxdepth 4 -type f -print || true

if [ -f "${JOB_OUTPUT}/logfile.log" ]; then
    cp "${JOB_OUTPUT}/logfile.log" "analysis_output/${TAG}/job_${JOBID}.pocket_coffea.log"
fi
if [ -d "${JOB_OUTPUT}/error" ]; then
    mkdir -p "analysis_output/${TAG}/job_${JOBID}.error"
    cp -r "${JOB_OUTPUT}/error/." "analysis_output/${TAG}/job_${JOBID}.error/"
fi
if [ -f "${JOB_OUTPUT}/failed_jobs.json" ]; then
    cp "${JOB_OUTPUT}/failed_jobs.json" "analysis_output/${TAG}/job_${JOBID}.failed_jobs.json"
fi

shopt -s nullglob
mapfile -t outputs < <(find "${JOB_OUTPUT}" -name '*.coffea' -type f -print)
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
