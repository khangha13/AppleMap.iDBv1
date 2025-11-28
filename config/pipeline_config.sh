#!/bin/bash
# =============================================================================
# GATK PIPELINE KH V1 - CENTRAL CONFIGURATION
# =============================================================================
# This file normalises configuration across the pipeline. Environment-specific
# values should live in config/environment.sh (generated from the template).
# The pipeline sources this module from every step.
# =============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [ -n "${PIPELINE_ROOT:-}" ]; then
    if [ -d "${PIPELINE_ROOT}" ]; then
        PIPELINE_ROOT="$(cd "${PIPELINE_ROOT}" && pwd)"
    else
        echo "[pipeline_config] ⚠️  Provided PIPELINE_ROOT (${PIPELINE_ROOT}) does not exist; falling back to config-relative path." >&2
        PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
    fi
else
    PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
fi
export PIPELINE_ROOT

if [ "${PIPELINE_CONFIG_SOURCED:-false}" = "true" ]; then
    return 0 2>/dev/null || exit 0
fi
PIPELINE_CONFIG_SOURCED="true"
export PIPELINE_CONFIG_SOURCED

ENV_FILE="${SCRIPT_DIR}/environment.sh"
ENV_TEMPLATE="${SCRIPT_DIR}/environment.template.sh"

# Ensure paths are absolute
ENV_FILE="$(cd "$(dirname "${ENV_FILE}")" && pwd)/$(basename "${ENV_FILE}")"
ENV_TEMPLATE="$(cd "$(dirname "${ENV_TEMPLATE}")" && pwd)/$(basename "${ENV_TEMPLATE}")"

if [ -f "${ENV_FILE}" ]; then
    source "${ENV_FILE}"
    PIPELINE_ENV_SOURCE="${ENV_FILE}"
elif [ -f "${ENV_TEMPLATE}" ]; then
    source "${ENV_TEMPLATE}"
    PIPELINE_ENV_SOURCE="${ENV_TEMPLATE}"
    echo "[pipeline_config] ⚠️  Using environment.template.sh; copy to environment.sh and update values for production." >&2
else
    echo "[pipeline_config] ❌ Missing configuration template at ${ENV_TEMPLATE}" >&2
    echo "[pipeline_config] ❌ SCRIPT_DIR: ${SCRIPT_DIR}" >&2
    exit 1
fi

# -----------------------------------------------------------------------------
# Helper utilities
# -----------------------------------------------------------------------------

_config_warn() {
    echo "[pipeline_config] ⚠️  $1" >&2
}

_config_info() {
    echo "[pipeline_config] ℹ️  $1" >&2
}

# Debug: Verify environment variables are loaded (only if verbose logging is enabled)
if [ "${PIPELINE_VERBOSE_LOGGING:-false}" = "true" ]; then
    _config_info "Environment source: ${PIPELINE_ENV_SOURCE}"
    _config_info "PIPELINE_REFERENCE_FASTA: ${PIPELINE_REFERENCE_FASTA:-<not set>}"
    _config_info "PIPELINE_ADAPTER_FASTA: ${PIPELINE_ADAPTER_FASTA:-<not set>}"
    _config_info "PIPELINE_KNOWN_SITES_VCF: ${PIPELINE_KNOWN_SITES_VCF:-<not set>}"
fi

_with_default() {
    local var_name="$1"
    local default_value="$2"
    local current_value="${!var_name:-}"
    if [ -z "${current_value}" ]; then
        printf -v "${var_name}" '%s' "${default_value}"
        _config_warn "${var_name} not set; defaulting to '${default_value}'"
    fi
}

# -----------------------------------------------------------------------------
# Metadata
# -----------------------------------------------------------------------------
PIPELINE_VERSION="1.0"
PIPELINE_NAME="GATK Pipeline KH v1"
PIPELINE_AUTHORS="Phu Khang Ha"
PIPELINE_CONTRIBUTORS="Paulo Henrique da Silva, Lisa Pasipanodya, Shashi Goonetilleke, Daniel Edge-Garza, Elizabeth Ross"
PIPELINE_DATE="October 2025"
PIPELINE_INSTITUTION="The University of Queensland"

# -----------------------------------------------------------------------------
# Storage paths
# -----------------------------------------------------------------------------
RDM_BASE_PATH="${PIPELINE_RDM_BASE:-${RDM_BASE_PATH:-/QRISdata/PLEASE_SET}}"
RDM_DATASETS_PATH="${PIPELINE_RDM_DATASETS:-${RDM_DATASETS_PATH:-${RDM_BASE_PATH}/WGS_Reference_Panel}}"
SCRATCH_BASE_PATH="${PIPELINE_SCRATCH_BASE:-${SCRATCH_BASE_PATH:-/scratch/user/$(whoami)}}"
LOG_BASE_PATH="${PIPELINE_LOG_BASE:-${LOG_BASE_PATH:-${SCRATCH_BASE_PATH}/logs}}"
WRAPPER_LOG_PATH="${WRAPPER_LOG_PATH:-${PIPELINE_ROOT}/logs/wrapper}"
PIPELINE_SLURM_SCRIPT_DIR="${PIPELINE_SLURM_SCRIPT_DIR:-${SCRATCH_BASE_PATH%/}/gatk_pipeline/slurm_scripts}"
PIPELINE_WORK_DIR="${PIPELINE_WORK_DIR:-${SCRATCH_BASE_PATH%/}/gatk_pipeline/work}"

# -----------------------------------------------------------------------------
# Reference & resource paths
# -----------------------------------------------------------------------------
# Initialize REFERENCE_GENOME safely (works with set -u)
# Always prioritize PIPELINE_REFERENCE_FASTA if set (handles multiple sourcing)
if [ -n "${PIPELINE_REFERENCE_FASTA:-}" ]; then
    REFERENCE_GENOME="${PIPELINE_REFERENCE_FASTA}"
    if [ "${PIPELINE_VERBOSE_LOGGING:-false}" = "true" ]; then
        _config_info "Set REFERENCE_GENOME from PIPELINE_REFERENCE_FASTA: ${REFERENCE_GENOME}"
    fi
elif [ -n "${REFERENCE_GENOME:-}" ] && [ "${REFERENCE_GENOME:-}" != "/path/to/reference.fasta" ]; then
    # REFERENCE_GENOME already set to a non-default value, keep it
    if [ "${PIPELINE_VERBOSE_LOGGING:-false}" = "true" ]; then
        _config_info "Using existing REFERENCE_GENOME: ${REFERENCE_GENOME}"
    fi
else
    REFERENCE_GENOME="/path/to/reference.fasta"
    if [ "${PIPELINE_VERBOSE_LOGGING:-false}" = "true" ]; then
        _config_warn "REFERENCE_GENOME not set, using default: ${REFERENCE_GENOME}"
    fi
fi

# If PIPELINE_REFERENCE_DIR is set and REFERENCE_GENOME is still the default,
# infer the path from the directory
if [[ -n "${PIPELINE_REFERENCE_DIR:-}" ]]; then
    if [ "${REFERENCE_GENOME:-}" = "/path/to/reference.fasta" ]; then
        REFERENCE_GENOME="${PIPELINE_REFERENCE_DIR%/}/reference.fasta"
        _config_warn "Inferring reference FASTA from PIPELINE_REFERENCE_DIR (${PIPELINE_REFERENCE_DIR})"
    fi
fi

# Validate REFERENCE_GENOME is set to a real path (warn if still default)
if [ "${REFERENCE_GENOME:-}" = "/path/to/reference.fasta" ]; then
    _config_warn "REFERENCE_GENOME is still set to default path. Check your environment configuration."
    _config_warn "Expected PIPELINE_REFERENCE_FASTA or PIPELINE_REFERENCE_DIR in ${PIPELINE_ENV_SOURCE}"
    _config_warn "PIPELINE_REFERENCE_FASTA value: ${PIPELINE_REFERENCE_FASTA:-<not set>}"
    _config_warn "PIPELINE_REFERENCE_DIR value: ${PIPELINE_REFERENCE_DIR:-<not set>}"
fi

REFERENCE_GENOME_INDEX="${REFERENCE_GENOME_INDEX:-${REFERENCE_GENOME:-/path/to/reference.fasta}.fai}"
# Use a safe reference for pattern removal
_ref_genome="${REFERENCE_GENOME:-/path/to/reference.fasta}"
REFERENCE_GENOME_DICT="${REFERENCE_GENOME_DICT:-${_ref_genome%.*}.dict}"
unset _ref_genome

KNOWN_SITES_VCF="${PIPELINE_KNOWN_SITES_VCF:-${KNOWN_SITES_VCF:-/path/to/known_sites.vcf.gz}}"
KNOWN_SITES_INDEX="${KNOWN_SITES_INDEX:-${KNOWN_SITES_VCF:-/path/to/known_sites.vcf.gz}.tbi}"

ADAPTER_FILE="${PIPELINE_ADAPTER_FASTA:-${ADAPTER_FILE:-/path/to/adapters/TruSeq2-PE.fa}}"

# -----------------------------------------------------------------------------
# SLURM defaults
# -----------------------------------------------------------------------------
PIPELINE_SLURM_ACCOUNT="${PIPELINE_SLURM_ACCOUNT:-a_qaafi_cas}"
PIPELINE_SLURM_PARTITION="${PIPELINE_SLURM_PARTITION:-general}"
PIPELINE_SLURM_QOS="${PIPELINE_SLURM_QOS:-}"
PIPELINE_SLURM_NODES="${PIPELINE_SLURM_NODES:-1}"
PIPELINE_SLURM_NTASKS="${PIPELINE_SLURM_NTASKS:-1}"

SLURM_ACCOUNT="${PIPELINE_SLURM_ACCOUNT}"
SLURM_PARTITION="${PIPELINE_SLURM_PARTITION}"
SLURM_QOS="${PIPELINE_SLURM_QOS}"
SLURM_NODES="${PIPELINE_SLURM_NODES}"
SLURM_NTASKS="${PIPELINE_SLURM_NTASKS}"

# -----------------------------------------------------------------------------
# Step 1A defaults
# -----------------------------------------------------------------------------
STEP1A_CPUS_PER_TASK="${STEP1A_CPUS:-10}"
STEP1A_MEMORY="${STEP1A_MEMORY:-32G}"
STEP1A_TIME_LIMIT="${STEP1A_TIME:-200:00:00}"
STEP1A_ARRAY_MAX="${STEP1A_ARRAY_LIMIT:-100}"
STEP1A_ACCOUNT="${STEP1A_ACCOUNT:-$PIPELINE_SLURM_ACCOUNT}"
STEP1A_PARTITION="${STEP1A_PARTITION:-$PIPELINE_SLURM_PARTITION}"
STEP1A_NODES="${STEP1A_NODES:-$PIPELINE_SLURM_NODES}"
STEP1A_NTASKS="${STEP1A_NTASKS:-$PIPELINE_SLURM_NTASKS}"

# -----------------------------------------------------------------------------
# Step 1B defaults
# -----------------------------------------------------------------------------
STEP1B_CPUS_PER_TASK="${STEP1B_CPUS:-6}"
STEP1B_MEMORY="${STEP1B_MEMORY:-36G}"
STEP1B_TIME_LIMIT="${STEP1B_TIME:-200:00:00}"
STEP1B_ARRAY_MAX="${STEP1B_ARRAY_LIMIT:-25}"
STEP1B_ACCOUNT="${STEP1B_ACCOUNT:-$PIPELINE_SLURM_ACCOUNT}"
STEP1B_PARTITION="${STEP1B_PARTITION:-$PIPELINE_SLURM_PARTITION}"
STEP1B_NODES="${STEP1B_NODES:-$PIPELINE_SLURM_NODES}"
STEP1B_NTASKS="${STEP1B_NTASKS:-$PIPELINE_SLURM_NTASKS}"

# -----------------------------------------------------------------------------
# Reference chromosome list (apple-specific default)
# -----------------------------------------------------------------------------
PIPELINE_CHROMOSOME_LIST="${PIPELINE_CHROMOSOME_LIST:-Chr00 Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11 Chr12 Chr13 Chr14 Chr15 Chr16 Chr17}"
# shellcheck disable=SC2206
PIPELINE_CHROMOSOMES=(${PIPELINE_CHROMOSOME_LIST})

# -----------------------------------------------------------------------------
# Step 1C defaults (placeholder for Beagle integration)
# -----------------------------------------------------------------------------
STEP1C_CPUS_PER_TASK="${STEP1C_CPUS:-8}"
STEP1C_MEMORY="${STEP1C_MEMORY:-48G}"
STEP1C_TIME_LIMIT="${STEP1C_TIME:-48:00:00}"
STEP1C_ACCOUNT="${STEP1C_ACCOUNT:-$PIPELINE_SLURM_ACCOUNT}"
STEP1C_PARTITION="${STEP1C_PARTITION:-$PIPELINE_SLURM_PARTITION}"
STEP1C_NODES="${STEP1C_NODES:-$PIPELINE_SLURM_NODES}"
STEP1C_NTASKS="${STEP1C_NTASKS:-$PIPELINE_SLURM_NTASKS}"

# -----------------------------------------------------------------------------
# Step 1D defaults
# -----------------------------------------------------------------------------
STEP1D_CPUS_PER_TASK="${STEP1D_CPUS:-16}"
STEP1D_MEMORY="${STEP1D_MEMORY:-256G}"
STEP1D_TIME_LIMIT="${STEP1D_TIME:-72:00:00}"
STEP1D_ACCOUNT="${STEP1D_ACCOUNT:-$PIPELINE_SLURM_ACCOUNT}"
STEP1D_PARTITION="${STEP1D_PARTITION:-$PIPELINE_SLURM_PARTITION}"
STEP1D_NODES="${STEP1D_NODES:-$PIPELINE_SLURM_NODES}"
STEP1D_NTASKS="${STEP1D_NTASKS:-$PIPELINE_SLURM_NTASKS}"
PLINK2_BIN="${PLINK2_BIN:-plink2}"
BCFTOOLS_BIN="${BCFTOOLS_BIN:-bcftools}"
STEP1D_PCA_ONLY="${STEP1D_PCA_ONLY:-false}"
STEP1D_REMOVE_RELATIVES="${STEP1D_REMOVE_RELATIVES:-false}"
STEP1D_PCA_DIR="${STEP1D_PCA_DIR:-pca_analysis}"
STEP1D_PCA_SHOW_LABELS="${STEP1D_PCA_SHOW_LABELS:-true}"
STEP1D_PCA_LABEL_SIZE="${STEP1D_PCA_LABEL_SIZE:-1.5}"
STEP1D_PCA_USE_GGREPEL="${STEP1D_PCA_USE_GGREPEL:-true}"

# -----------------------------------------------------------------------------
# Software paths (allow module system to populate)
# -----------------------------------------------------------------------------
GATK_COMMAND="${GATK_COMMAND:-gatk}"
if [ -z "${TRIMMOMATIC_JAR:-}" ]; then
    if [ -n "${EBROOTTRIMMOMATIC:-}" ]; then
        TRIMMOMATIC_JAR="${EBROOTTRIMMOMATIC}/trimmomatic-0.39.jar"
    else
        TRIMMOMATIC_JAR="trimmomatic-0.39.jar"
    fi
fi
BWA_PATH="${BWA_PATH:-bwa}"
SAMTOOLS_PATH="${SAMTOOLS_PATH:-samtools}"
FASTQC_PATH="${FASTQC_PATH:-fastqc}"

# -----------------------------------------------------------------------------
# Pipeline parameters
# -----------------------------------------------------------------------------
# Trimmomatic parameters (matching original Apple_GATK_Pipeline_1a_CallVariantsPerSampleKHv2.sh)
TRIMMOMATIC_LEADING="${TRIMMOMATIC_LEADING:-20}"
TRIMMOMATIC_TRAILING="${TRIMMOMATIC_TRAILING:-20}"
TRIMMOMATIC_SLIDINGWINDOW="${TRIMMOMATIC_SLIDINGWINDOW:-3:15}"
TRIMMOMATIC_AVGQUAL="${TRIMMOMATIC_AVGQUAL:-20}"
TRIMMOMATIC_MINLEN="${TRIMMOMATIC_MINLEN:-36}"
TRIMMOMATIC_ILLUMINACLIP_SETTINGS="${TRIMMOMATIC_ILLUMINACLIP_SETTINGS:-2:30:3:1:True}"

BWA_MEM_FLAGS="${BWA_MEM_FLAGS:- -M}"
GATK_BATCH_SIZE="${GATK_BATCH_SIZE:-50}"
GATK_READER_THREADS="${GATK_READER_THREADS:-4}"

# -----------------------------------------------------------------------------
# Logging configuration
# -----------------------------------------------------------------------------
DEFAULT_LOG_LEVEL="${DEFAULT_LOG_LEVEL:-INFO}"
PIPELINE_LOG_LEVEL="${PIPELINE_LOG_LEVEL:-${DEFAULT_LOG_LEVEL}}"
LOG_TO_CONSOLE="${LOG_TO_CONSOLE:-true}"
LOG_TO_FILE="${LOG_TO_FILE:-true}"
PIPELINE_VERBOSE_LOGGING="${PIPELINE_VERBOSE_LOGGING:-false}"
PIPELINE_GLOBAL_DRY_RUN="${PIPELINE_GLOBAL_DRY_RUN:-false}"

# -----------------------------------------------------------------------------
# Backup and shared reference configuration
# -----------------------------------------------------------------------------
PIPELINE_ENABLE_BACKUP="${PIPELINE_ENABLE_BACKUP:-true}"
PIPELINE_SHARED_REF_TIMEOUT="${PIPELINE_SHARED_REF_TIMEOUT:-120}"
PIPELINE_SHARED_REF_POLL_INTERVAL="${PIPELINE_SHARED_REF_POLL_INTERVAL:-3}"
PIPELINE_BACKUP_STEPS="${PIPELINE_BACKUP_STEPS:-3,4,5,6}"

# -----------------------------------------------------------------------------
# Validation configuration
# -----------------------------------------------------------------------------
REQUIRED_DIRS=(
    "1.FASTQ"
    "2.FASTQC_pre_trimmed"
    "3.FASTQC_post_trimmed"
    "4.BAM"
    "5.Individual_VCF"
    "6.genomicsdb_output"
    "7.Consolidated_VCF"
    "8.Imputated_VCF_BEAGLE"
    "9.Imputation_QC"
)

FASTQ_EXTENSIONS=("_1.fastq.gz" "_2.fastq.gz")
GVCF_EXTENSIONS=("_raw.g.vcf.gz" "_raw.g.vcf.gz.tbi")
VCF_EXTENSIONS=("_genotyped.vcf.gz" "_genotyped.vcf.gz.tbi")
CONSOLIDATED_VCF_EXTENSIONS=(".vcf.gz" ".vcf.gz.tbi")

# -----------------------------------------------------------------------------
# Monitoring configuration
# -----------------------------------------------------------------------------
MONITOR_INTERVAL="${MONITOR_INTERVAL:-30}"
MAX_MONITOR_ATTEMPTS="${MAX_MONITOR_ATTEMPTS:-100}"
EMAIL_NOTIFICATIONS="${EMAIL_NOTIFICATIONS:-false}"
EMAIL_ADDRESS="${EMAIL_ADDRESS:-}"

# -----------------------------------------------------------------------------
# Advanced configuration
# -----------------------------------------------------------------------------
ENABLE_RESUME="${ENABLE_RESUME:-true}"
RESUME_STEP_FILE="${RESUME_STEP_FILE:-resume_step.txt}"
ENABLE_BACKUPS="${ENABLE_BACKUPS:-true}"
BACKUP_INTERVAL="${BACKUP_INTERVAL:-5}"
CLEANUP_TEMP_FILES="${CLEANUP_TEMP_FILES:-true}"
CLEANUP_INTERMEDIATE_FILES="${CLEANUP_INTERMEDIATE_FILES:-false}"
DEBUG_MODE="${DEBUG_MODE:-false}"
TEST_MODE="${TEST_MODE:-false}"
VERBOSE_OUTPUT="${VERBOSE_OUTPUT:-false}"

# -----------------------------------------------------------------------------
# Helper accessors
# -----------------------------------------------------------------------------
get_step1a_config() {
    echo "CPUS=${STEP1A_CPUS_PER_TASK}"
    echo "MEMORY=${STEP1A_MEMORY}"
    echo "TIME=${STEP1A_TIME_LIMIT}"
    echo "ACCOUNT=${STEP1A_ACCOUNT}"
    echo "PARTITION=${STEP1A_PARTITION}"
    echo "NODES=${STEP1A_NODES}"
    echo "NTASKS=${STEP1A_NTASKS}"
    echo "ARRAY_MAX=${STEP1A_ARRAY_MAX}"
    [ -n "${SLURM_QOS}" ] && echo "QOS=${SLURM_QOS}"
}

get_step1b_config() {
    echo "CPUS=${STEP1B_CPUS_PER_TASK}"
    echo "MEMORY=${STEP1B_MEMORY}"
    echo "TIME=${STEP1B_TIME_LIMIT}"
    echo "ACCOUNT=${STEP1B_ACCOUNT}"
    echo "PARTITION=${STEP1B_PARTITION}"
    echo "NODES=${STEP1B_NODES}"
    echo "NTASKS=${STEP1B_NTASKS}"
    echo "ARRAY_MAX=${STEP1B_ARRAY_MAX}"
    [ -n "${SLURM_QOS}" ] && echo "QOS=${SLURM_QOS}"
}

get_step1c_config() {
    echo "CPUS=${STEP1C_CPUS_PER_TASK}"
    echo "MEMORY=${STEP1C_MEMORY}"
    echo "TIME=${STEP1C_TIME_LIMIT}"
    echo "ACCOUNT=${STEP1C_ACCOUNT}"
    echo "PARTITION=${STEP1C_PARTITION}"
    echo "NODES=${STEP1C_NODES}"
    echo "NTASKS=${STEP1C_NTASKS}"
    [ -n "${SLURM_QOS}" ] && echo "QOS=${SLURM_QOS}"
}

get_step1d_config() {
    echo "CPUS=${STEP1D_CPUS_PER_TASK}"
    echo "MEMORY=${STEP1D_MEMORY}"
    echo "TIME=${STEP1D_TIME_LIMIT}"
    echo "ACCOUNT=${STEP1D_ACCOUNT}"
    echo "PARTITION=${STEP1D_PARTITION}"
    echo "NODES=${STEP1D_NODES}"
    echo "NTASKS=${STEP1D_NTASKS}"
    [ -n "${SLURM_QOS}" ] && echo "QOS=${SLURM_QOS}"
}

get_reference_fasta() {
    echo "${REFERENCE_GENOME:-/path/to/reference.fasta}"
}

get_known_sites_vcf() {
    echo "${KNOWN_SITES_VCF:-/path/to/known_sites.vcf.gz}"
}

get_adapter_fasta() {
    echo "${ADAPTER_FILE:-/path/to/adapters/TruSeq2-PE.fa}"
}

get_rdm_datasets_path() {
    echo "${RDM_DATASETS_PATH}"
}

# Legacy alias for backward compatibility (deprecated)
get_rdm_base_path() {
    get_rdm_datasets_path
}

get_scratch_base_path() {
    echo "${SCRATCH_BASE_PATH}"
}

get_log_base_path() {
    echo "${LOG_BASE_PATH}"
}

get_pipeline_env_source() {
    echo "${PIPELINE_ENV_SOURCE}"
}

# Provide a one-line summary for debugging
pipeline_config_summary() {
    _config_info "Environment source: ${PIPELINE_ENV_SOURCE}"
    _config_info "RDM base path: ${RDM_BASE_PATH}"
    _config_info "RDM datasets path: ${RDM_DATASETS_PATH}"
    _config_info "Scratch base path: ${SCRATCH_BASE_PATH}"
    _config_info "Reference FASTA: ${REFERENCE_GENOME:-/path/to/reference.fasta}"
}
