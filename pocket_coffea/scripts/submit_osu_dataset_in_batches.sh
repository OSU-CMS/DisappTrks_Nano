#!/bin/bash
set -euo pipefail

if [ "$#" -lt 4 ]; then
    cat <<'USAGE' >&2
Usage:
  scripts/submit_osu_dataset_in_batches.sh DATASET_JSON TAG SAMPLE YEAR

This is a convenience wrapper around submit_osu_dataset.sh for the OSU mode
where worker nodes cannot see the ROOT-file storage path.  It submits one ROOT
file per job with Condor file transfer enabled, but groups submissions into
offset batches so a transient condor_submit hiccup does not ruin the day.

Common environment overrides:
  BATCH_SIZE       default: 100
  START_OFFSET     default: 0
  STOP_OFFSET      default: all jobs
  SUBMIT_PAUSE     default: 30; seconds between batches
  CATEGORY_MODE    default: muon_pveto
  CHUNKSIZE        default: 50000
  JOB_TIMEOUT      default: 90m
  REQUEST_MEMORY   default: 4GB
  REQUEST_DISK     default: 8GB
  SUBMIT_RETRIES   default: 10
  SUBMIT_SLEEP     default: 1

Example:
  CATEGORY_MODE=electron_pveto BATCH_SIZE=100 \
    scripts/submit_osu_dataset_in_batches.sh \
      datasets/osu_2022CD_EGamma.json 2022CD_electron_pveto_v2 DATA_EGamma 2022_preEE
USAGE
    exit 2
fi

DATASET_JSON="$1"
TAG="$2"
SAMPLE="$3"
YEAR="$4"

BATCH_SIZE="${BATCH_SIZE:-100}"
START_OFFSET="${START_OFFSET:-0}"
STOP_OFFSET="${STOP_OFFSET:-}"
SUBMIT_PAUSE="${SUBMIT_PAUSE:-30}"

if [ "${BATCH_SIZE}" -le 0 ]; then
    echo "BATCH_SIZE must be positive: ${BATCH_SIZE}" >&2
    exit 1
fi
if [ "${START_OFFSET}" -lt 0 ]; then
    echo "START_OFFSET must be non-negative: ${START_OFFSET}" >&2
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

FILES_PER_JOB="${FILES_PER_JOB:-1}"
if [ "${FILES_PER_JOB}" -ne 1 ]; then
    echo "This wrapper requires FILES_PER_JOB=1 for ROOT-file transfer mode." >&2
    exit 1
fi

TOTAL_JOBS="${NFILES}"
if [ -z "${STOP_OFFSET}" ] || [ "${STOP_OFFSET}" -gt "${TOTAL_JOBS}" ]; then
    STOP_OFFSET="${TOTAL_JOBS}"
fi
if [ "${STOP_OFFSET}" -le "${START_OFFSET}" ]; then
    echo "Nothing to submit: START_OFFSET=${START_OFFSET}, STOP_OFFSET=${STOP_OFFSET}" >&2
    exit 1
fi

echo "Submitting ${TAG} in batches"
echo "  dataset:      ${DATASET_JSON}"
echo "  sample/year:  ${SAMPLE} / ${YEAR}"
echo "  total files:  ${NFILES}"
echo "  offsets:      ${START_OFFSET}..$((STOP_OFFSET - 1))"
echo "  batch size:   ${BATCH_SIZE}"
echo "  pause:        ${SUBMIT_PAUSE}s"
echo "  category:     ${CATEGORY_MODE:-muon_pveto}"

offset="${START_OFFSET}"
while [ "${offset}" -lt "${STOP_OFFSET}" ]; do
    remaining=$(( STOP_OFFSET - offset ))
    this_batch="${BATCH_SIZE}"
    if [ "${remaining}" -lt "${this_batch}" ]; then
        this_batch="${remaining}"
    fi

    echo
    echo "Submitting batch offset=${offset}, jobs=${this_batch}"
    TRANSFER_ROOT_FILES=1 \
    FILES_PER_JOB=1 \
    MAX_JOBS="${this_batch}" \
    JOB_OFFSET="${offset}" \
    CHUNKSIZE="${CHUNKSIZE:-50000}" \
    JOB_TIMEOUT="${JOB_TIMEOUT:-90m}" \
    REQUEST_MEMORY="${REQUEST_MEMORY:-4GB}" \
    REQUEST_DISK="${REQUEST_DISK:-8GB}" \
    CATEGORY_MODE="${CATEGORY_MODE:-muon_pveto}" \
    SUBMIT_RETRIES="${SUBMIT_RETRIES:-10}" \
    SUBMIT_SLEEP="${SUBMIT_SLEEP:-1}" \
        scripts/submit_osu_dataset.sh "${DATASET_JSON}" "${TAG}" "${SAMPLE}" "${YEAR}"

    offset=$(( offset + this_batch ))
    if [ "${offset}" -lt "${STOP_OFFSET}" ] && [ "${SUBMIT_PAUSE}" != "0" ]; then
        echo "Sleeping ${SUBMIT_PAUSE}s before next batch..."
        sleep "${SUBMIT_PAUSE}"
    fi
done

echo
echo "Submitted requested batches for ${TAG}."
