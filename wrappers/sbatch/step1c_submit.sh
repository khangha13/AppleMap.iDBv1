#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [ -n "${PIPELINE_ROOT:-}" ]; then
    if [ -d "${PIPELINE_ROOT}" ]; then
        PIPELINE_ROOT="$(cd "${PIPELINE_ROOT}" && pwd)"
    else
        echo "[step1c_submit] ⚠️  Provided PIPELINE_ROOT (${PIPELINE_ROOT}) does not exist; falling back to wrapper-relative path." >&2
        PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
    fi
else
    PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
fi
export PIPELINE_ROOT

usage() {
    cat <<EOF
Usage: ${0##*/} <dataset_name> [rdm_base_path]

Submits the Step 1C Beagle imputation to Slurm with resources derived from config/pipeline_config.sh.

Arguments:
  dataset_name    Dataset identifier (required)
  rdm_base_path   RDM dataset path (optional - will be inferred from config if omitted)

If rdm_base_path is omitted, it will be automatically inferred as:
  \${RDM_DATASETS_PATH}/\${dataset_name}

Examples:
  ${0##*/} MyDataset
  ${0##*/} MyDataset /QRISdata/Q8367/WGS_Reference_Panel/MyDataset
EOF
    exit 1
}

if [ $# -lt 1 ]; then
    usage
fi

DATASET_NAME="$1"

# Guard: reject flag-like strings consumed as positional arguments.
if [[ "${DATASET_NAME}" == -* ]]; then
    echo "[step1c_submit] ❌ First argument '${DATASET_NAME}' looks like a flag, not a dataset name." >&2
    echo "[step1c_submit] Expected usage: ${0##*/} <dataset_name> [rdm_base_path]" >&2
    exit 1
fi

# Source config early to get helper functions
source "${PIPELINE_ROOT}/config/pipeline_config.sh"

# Infer RDM path if not provided (similar to interactive wrappers)
if [ $# -eq 1 ]; then
    DATASET_RDM_PATH="$(get_rdm_datasets_path)/${DATASET_NAME}"
    echo "[step1c_submit] ℹ️  Auto-inferred RDM path: ${DATASET_RDM_PATH}"
elif [ $# -eq 2 ]; then
    DATASET_RDM_PATH="$2"
    echo "[step1c_submit] Using provided RDM path: ${DATASET_RDM_PATH}"
else
    usage
fi

# Pass to the module runner (similar to step1a/step1d pattern)
exec bash "${PIPELINE_ROOT}/modules/step1c/bin/run_step1c.sh" "${DATASET_NAME}" "${DATASET_RDM_PATH}"
