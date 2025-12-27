#!/bin/bash
# STEP 1A ARRAY WORKER TEMPLATE
set -euo pipefail

DATASET_PATH="$1"
SAMPLE_LIST_FILE="$2"

if [ -z "${DATASET_PATH:-}" ] || [ -z "${SAMPLE_LIST_FILE:-}" ]; then
    echo "Usage: ${BASH_SOURCE[0]} <dataset_path> <sample_list_file>" >&2
    exit 1
fi

# =============================================================================
# MODULE LOADING SECTION
# =============================================================================
# Load all required bioinformatics tools and their specific versions
# These versions have been tested and are compatible with each other

module purge
module load fastqc/0.11.9-java-11              # Quality control for raw reads
module load trimmomatic/0.39-java-11           # Read trimming and adapter removal
module load bwa/0.7.17-gcccore-11.3.0          # Burrows-Wheeler Aligner for read mapping
module load samtools/1.16.1-gcc-11.3.0         # SAM/BAM file manipulation and indexing
module load gatk/4.3.0.0-gcccore-11.3.0-java-11 # Genome Analysis Toolkit for variant calling
module load picard                              # Java tools for working with sequencing data

if [ -n "${PIPELINE_ROOT:-}" ]; then
    if [ -d "${PIPELINE_ROOT}" ]; then
        PIPELINE_ROOT="$(cd "${PIPELINE_ROOT}" && pwd)"
    else
        echo "[step1a_array] ⚠️  Provided PIPELINE_ROOT (${PIPELINE_ROOT}) does not exist; falling back to template-relative path." >&2
        PIPELINE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
    fi
else
    PIPELINE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
fi
export PIPELINE_ROOT
STEP1A_MODULE_DIR="${PIPELINE_ROOT}/modules/step1a"

source "${PIPELINE_ROOT}/config/pipeline_config.sh"
source "${PIPELINE_ROOT}/lib/pipeline_common.sh"
source "${STEP1A_MODULE_DIR}/lib/functions.sh"
source "${STEP1A_MODULE_DIR}/bin/run_step1a.sh"

# Initialize logging for the array job
init_logging "step1a" "pipeline" "${DATASET_NAME}"

if [ ! -f "${SAMPLE_LIST_FILE}" ]; then
    error_exit "Sample list file not found: ${SAMPLE_LIST_FILE}"
fi

mapfile -t STEP1A_SAMPLES < "${SAMPLE_LIST_FILE}"

if [ -z "${SLURM_ARRAY_TASK_ID:-}" ]; then
    error_exit "SLURM_ARRAY_TASK_ID is not set; this script must run as an array job."
fi
n_samples="${#STEP1A_SAMPLES[@]}"
if ! [[ "${SLURM_ARRAY_TASK_ID}" =~ ^[0-9]+$ ]] || [ "${SLURM_ARRAY_TASK_ID}" -lt 0 ] || [ "${SLURM_ARRAY_TASK_ID}" -ge "${n_samples}" ]; then
    error_exit "Array index ${SLURM_ARRAY_TASK_ID} out of range 0..$((n_samples-1))"
fi

SAMPLE="${STEP1A_SAMPLES[$SLURM_ARRAY_TASK_ID]}"

if [ -z "${SAMPLE:-}" ]; then
    error_exit "No sample found for array index ${SLURM_ARRAY_TASK_ID}"
fi

log_info "Step 1A array task ${SLURM_ARRAY_TASK_ID} processing sample ${SAMPLE}"

# Extract dataset name from DATASET_PATH for shared reference setup
DATASET_NAME="$(basename "${DATASET_PATH}")"

ref_genome="$(get_reference_fasta)"
known_sites="$(get_known_sites_vcf)"
adapter_file="$(get_adapter_fasta)"

ensure_shared_references_ready "${DATASET_NAME}" "${ref_genome}" "${known_sites}" "${adapter_file}"

# Task 0 can still initialize shared backup directories
if [ "${SLURM_ARRAY_TASK_ID}" -eq 0 ]; then
    log_info "Task 0: Initializing backup directories"
    backup_base="${SCRATCH_BASE_PATH%/}/${DATASET_NAME}_backup"
    mkdir -p "${backup_base}"
    log_info "Backup directories initialized"
fi

execute_step1a_pipeline "${SAMPLE}" "${DATASET_PATH}" "${DATASET_NAME}"
