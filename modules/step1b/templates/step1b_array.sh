#!/bin/bash
# STEP 1B ARRAY WORKER TEMPLATE
set -euo pipefail

DATASET_PATH="$1"
CHROMOSOME_LIST_FILE="$2"

if [ -z "${DATASET_PATH:-}" ] || [ -z "${CHROMOSOME_LIST_FILE:-}" ]; then
    echo "Usage: ${BASH_SOURCE[0]} <dataset_path> <chromosome_list_file>" >&2
    exit 1
fi

step1b_array_exit_trap() {
    local status=$?
    if [ $status -ne 0 ] && [ -n "${STEP1B_FAILURE_FLAG_PATH:-}" ]; then
        if [ ! -f "${STEP1B_FAILURE_FLAG_PATH}" ]; then
            {
                printf '%s\tStep 1B array task %s (chromosome %s) failed with status %s\n' \
                    "$(date +%Y-%m-%dT%H:%M:%S)" \
                    "${SLURM_ARRAY_TASK_ID:-unknown}" \
                    "${chromosome:-unknown}" \
                    "$status"
            } > "${STEP1B_FAILURE_FLAG_PATH}" 2>/dev/null || true
        fi
    fi
    exit $status
}

# =============================================================================
# MODULE LOADING SECTION
# =============================================================================
# Load all required bioinformatics tools and their specific versions
# These versions have been tested and are compatible with each other

module purge
module load gatk/4.3.0.0-gcccore-11.3.0-java-11 # Genome Analysis Toolkit for variant calling
module load samtools/1.16.1-gcc-11.3.0         # SAM/BAM file manipulation and indexing

PIPELINE_ROOT="__PIPELINE_ROOT_PLACEHOLDER__"
export PIPELINE_ROOT
STEP1B_MODULE_DIR="${PIPELINE_ROOT}/modules/step1b"

source "${PIPELINE_ROOT}/config/pipeline_config.sh"
source "${PIPELINE_ROOT}/lib/pipeline_common.sh"
source "${STEP1B_MODULE_DIR}/lib/functions.sh"
source "${STEP1B_MODULE_DIR}/bin/run_step1b.sh"

# Initialize logging for the array job
init_logging "step1b" "pipeline" "${DATASET_NAME}"

if [ ! -f "${CHROMOSOME_LIST_FILE}" ]; then
    error_exit "Chromosome list file not found: ${CHROMOSOME_LIST_FILE}"
fi

mapfile -t STEP1B_CHROMS < "${CHROMOSOME_LIST_FILE}"

if [ -z "${SLURM_ARRAY_TASK_ID:-}" ]; then
    error_exit "SLURM_ARRAY_TASK_ID is not set; this script must run as an array job."
fi
n_chroms="${#STEP1B_CHROMS[@]}"
if ! [[ "${SLURM_ARRAY_TASK_ID}" =~ ^[0-9]+$ ]] || [ "${SLURM_ARRAY_TASK_ID}" -lt 0 ] || [ "${SLURM_ARRAY_TASK_ID}" -ge "${n_chroms}" ]; then
    error_exit "Array index ${SLURM_ARRAY_TASK_ID} out of range 0..$((n_chroms-1))"
fi

chromosome="${STEP1B_CHROMS[$SLURM_ARRAY_TASK_ID]}"
if [ -z "${chromosome:-}" ]; then
    error_exit "No chromosome found for array index ${SLURM_ARRAY_TASK_ID}"
fi

log_info "Step 1B array task ${SLURM_ARRAY_TASK_ID} processing ${chromosome}"
trap 'step1b_array_exit_trap' EXIT

# Extract dataset name from DATASET_PATH for shared reference setup
DATASET_NAME="$(basename "${DATASET_PATH}")"

ref_genome="$(get_reference_fasta)"

ensure_shared_references_ready "${DATASET_NAME}" "${ref_genome}"

execute_step1b_pipeline "${chromosome}" "${DATASET_PATH}" "${DATASET_NAME}"
