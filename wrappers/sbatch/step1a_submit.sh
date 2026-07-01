#!/bin/bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [ -n "${PIPELINE_ROOT:-}" ]; then
    if [ -d "${PIPELINE_ROOT}" ]; then
        PIPELINE_ROOT="$(cd "${PIPELINE_ROOT}" && pwd)"
    else
        echo "[step1a_submit] ⚠️  Provided PIPELINE_ROOT (${PIPELINE_ROOT}) does not exist; falling back to wrapper-relative path." >&2
        PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
    fi
else
    PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
fi
export PIPELINE_ROOT

usage() {
    cat <<EOF
Usage: ${0##*/} <dataset_name> [rdm_base_path] [step1a options]

Examples:
  ${0##*/} NCBI_truth_set --from-recal-bam --bam-pattern 'remainder_*.bam'
  ${0##*/} NCBI_truth_set /QRISdata/Q8367/WGS_Reference_Panel/NCBI_truth_set --from-recal-bam --bam-pattern 'remainder_*.bam'
EOF
    exit 1
}

if [ $# -lt 1 ]; then
    usage
fi

DATASET_NAME="$1"
shift

if [[ "${DATASET_NAME}" == -* ]]; then
    echo "[step1a_submit] ❌ First argument '${DATASET_NAME}' looks like a flag, not a dataset name." >&2
    usage
fi

source "${PIPELINE_ROOT}/config/pipeline_config.sh"

if [ $# -gt 0 ] && [[ "${1}" != -* ]]; then
    DATASET_RDM_PATH="$1"
    shift
    echo "[step1a_submit] Using provided RDM path: ${DATASET_RDM_PATH}"
else
    DATASET_RDM_PATH="$(get_rdm_datasets_path)/${DATASET_NAME}"
    echo "[step1a_submit] ℹ️  Auto-inferred RDM path: ${DATASET_RDM_PATH}"
fi

exec bash "${PIPELINE_ROOT}/modules/step1a/bin/run_step1a.sh" "${DATASET_NAME}" "${DATASET_RDM_PATH}" "$@"
