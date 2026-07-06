#!/usr/bin/env bash
# Short iterative tau-pveto smoke tests for 2022/2023 Muon and EGamma data.
#
# Run from DisappTrks_Nano/pocket_coffea:
#   scripts/run_tau_pveto_smokes.sh
#
# Useful overrides:
#   LIMIT_FILES=2 LIMIT_CHUNKS=2 scripts/run_tau_pveto_smokes.sh
#   CASES="2023C_mu 2023D_eg" scripts/run_tau_pveto_smokes.sh

set -euo pipefail

LIMIT_FILES="${LIMIT_FILES:-1}"
LIMIT_CHUNKS="${LIMIT_CHUNKS:-1}"
CHUNKSIZE="${CHUNKSIZE:-50000}"
OUTPUT_BASE="${OUTPUT_BASE:-analysis_output/tau_pveto_smoke}"
CASES="${CASES:-2022CD_mu 2022EFG_mu 2023C_mu 2023D_mu 2022CD_eg 2022EFG_eg 2023C_eg 2023D_eg}"

if [ ! -f "config.py" ] || [ ! -d "datasets" ]; then
    echo "ERROR: run this from DisappTrks_Nano/pocket_coffea" >&2
    exit 2
fi

# For local diagnostics, do not require local JME payloads. Disable this for
# production-like checks once payload access is confirmed.
export DISAPPTRKS_ALLOW_MISSING_JET_VETO_MAP="${DISAPPTRKS_ALLOW_MISSING_JET_VETO_MAP:-1}"
unset DISAPPTRKS_ENABLE_SEARCH_DIAGNOSTICS

run_case() {
    local label="$1"
    local dataset_json="$2"
    local category_mode="$3"
    local outputdir="${OUTPUT_BASE}/${label}"

    echo
    echo "==> ${label}: ${dataset_json}, ${category_mode}"
    DISAPPTRKS_DATASET_JSON="${dataset_json}" \
    DISAPPTRKS_CATEGORY_MODE="${category_mode}" \
    python -m pocket_coffea.scripts.runner run \
        --cfg config.py \
        --outputdir "${outputdir}" \
        --executor iterative \
        --limit-files "${LIMIT_FILES}" \
        --limit-chunks "${LIMIT_CHUNKS}" \
        --chunksize "${CHUNKSIZE}" \
        --process-separately
}

for case_name in ${CASES}; do
    case "${case_name}" in
        2022CD_mu)  run_case "${case_name}_tau_mu"  "datasets/eos_2022CD_Muon.json"   "tau_mu_pveto" ;;
        2022EFG_mu) run_case "${case_name}_tau_mu"  "datasets/eos_2022EFG_Muon.json"  "tau_mu_pveto" ;;
        2023C_mu)   run_case "${case_name}_tau_mu"  "datasets/eos_2023C_muon.json"    "tau_mu_pveto" ;;
        2023D_mu)   run_case "${case_name}_tau_mu"  "datasets/eos_2023D_Muon.json"    "tau_mu_pveto" ;;
        2022CD_eg)  run_case "${case_name}_tau_ele" "datasets/eos_2022CD_EGamma.json" "tau_ele_pveto" ;;
        2022EFG_eg) run_case "${case_name}_tau_ele" "datasets/eos_2022EFG_EGamma.json" "tau_ele_pveto" ;;
        2023C_eg)   run_case "${case_name}_tau_ele" "datasets/eos_2023C_EGamma.json"  "tau_ele_pveto" ;;
        2023D_eg)   run_case "${case_name}_tau_ele" "datasets/eos_2023D_EGamma.json"  "tau_ele_pveto" ;;
        *)
            echo "ERROR: unknown case '${case_name}'" >&2
            exit 2
            ;;
    esac
done
