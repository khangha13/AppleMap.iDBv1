#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [ -n "${PIPELINE_ROOT:-}" ]; then
    if [ -d "${PIPELINE_ROOT}" ]; then
        PIPELINE_ROOT="$(cd "${PIPELINE_ROOT}" && pwd)"
    else
        echo "[step1b_submit] ⚠️  Provided PIPELINE_ROOT (${PIPELINE_ROOT}) does not exist; falling back to wrapper-relative path." >&2
        PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
    fi
else
    PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
fi
export PIPELINE_ROOT

usage() {
    cat <<EOF
Usage: ${0##*/} <dataset_name> <rdm_base_path>

Submits the Step 1B orchestrator to Slurm with resources derived from config/pipeline_config.sh.
EOF
    exit 1
}

if [ $# -lt 2 ]; then
    usage
fi

DATASET_NAME="$1"
DATASET_RDM_PATH="$2"

shift 2 || true

source "${PIPELINE_ROOT}/config/pipeline_config.sh"

LOG_DIR="${MASTER_LOG_DIR:-${LOG_BASE_PATH%/}/${DATASET_NAME}}"
mkdir -p "${LOG_DIR}"

RUN_STAMP="$(date +%Y%m%d_%H%M%S)"
STDOUT_PATH="${LOG_DIR}/step1b_master_${RUN_STAMP}_%j.output"
STDERR_PATH="${LOG_DIR}/step1b_master_${RUN_STAMP}_%j.error"

SCRIPT_PATH="${PIPELINE_ROOT}/modules/step1b/bin/run_step1b.sh"
if [ ! -x "${SCRIPT_PATH}" ]; then
    echo "[step1b_submit] ❌ Step 1B runner not found or not executable: ${SCRIPT_PATH}" >&2
    exit 1
fi

JOB_NAME="GATK_step1b_${DATASET_NAME}"

sbatch_cmd=(
    sbatch
    --parsable
    -J "${JOB_NAME}"
    -A "${STEP1B_ACCOUNT}"
    -p "${STEP1B_PARTITION}"
    -N "${STEP1B_NODES}"
    -n "${STEP1B_NTASKS}"
    -c "${STEP1B_CPUS_PER_TASK}"
    --mem="${STEP1B_MEMORY}"
    -t "${STEP1B_TIME_LIMIT}"
    -o "${STDOUT_PATH}"
    -e "${STDERR_PATH}"
    --export=ALL,PIPELINE_ROOT="${PIPELINE_ROOT}"
)

if [ -n "${PIPELINE_SLURM_QOS:-}" ]; then
    sbatch_cmd+=( --qos "${PIPELINE_SLURM_QOS}" )
fi

echo "[step1b_submit] Dataset: ${DATASET_NAME}"
echo "[step1b_submit] RDM dataset path: ${DATASET_RDM_PATH}"
echo "[step1b_submit] Log directory: ${LOG_DIR}"

if ! sbatch_output=$("${sbatch_cmd[@]}" "${SCRIPT_PATH}" "${DATASET_NAME}" "${DATASET_RDM_PATH}"); then
    echo "[step1b_submit] ❌ Failed to submit Step 1B job via sbatch." >&2
    exit 1
fi

job_id="$(printf '%s\n' "${sbatch_output}" | head -n1 | tr -d '\r')"

if [ -z "${job_id}" ]; then
    echo "[step1b_submit] ⚠️  sbatch returned empty job ID (raw output: ${sbatch_output})." >&2
else
    echo "[step1b_submit] ✅ Step 1B job submitted (Job ID: ${job_id})."
    echo "[step1b_submit]    stdout: ${STDOUT_PATH//%j/${job_id}}"
    echo "[step1b_submit]    stderr: ${STDERR_PATH//%j/${job_id}}"
fi
