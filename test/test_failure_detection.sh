#!/bin/bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export PIPELINE_ROOT="${ROOT_DIR}"
export SCRATCH_BASE_PATH="${SCRATCH_BASE_PATH:-/tmp/gatk_scratch}"
mkdir -p "${SCRATCH_BASE_PATH}"

TMP_BASE="${TMPDIR:-/tmp}/gatk_failure_test"
mkdir -p "${TMP_BASE}"
export PIPELINE_LOG_BASE="${TMP_BASE}/logs"

source "${PIPELINE_ROOT}/config/pipeline_config.sh"
source "${PIPELINE_ROOT}/lib/logging.sh"
init_logging "test_failure_detection" "pipeline"
source "${PIPELINE_ROOT}/lib/validation.sh"

DATASET_DIR="$(mktemp -d "${TMP_BASE}/dataset_XXXX")"
mkdir -p "${DATASET_DIR}/7.Consolidated_VCF"
DATASET_NAME="$(basename "${DATASET_DIR}")"
FLAG_PATH="${PIPELINE_LOG_BASE}/${DATASET_NAME}/step1b_failed.flag"
mkdir -p "$(dirname "${FLAG_PATH}")"
echo "synthetic failure for test" > "${FLAG_PATH}"

status="$(check_step1b_status "${DATASET_DIR}")"
if [ "${status}" != "Failed" ]; then
    echo "❌ Expected Step 1B status 'Failed', got '${status}'"
    exit 1
fi

echo "✅ Step 1B failure flag detection behaves as expected."

