#!/bin/bash
# =============================================================================
# STEP 1B MODULE - COMBINE GVCFs
# =============================================================================

STEP1B_BIN_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Resolve PIPELINE_ROOT even when the script is executed from a SLURM spool copy
if [ -n "${PIPELINE_ROOT:-}" ] && [ -d "${PIPELINE_ROOT}" ]; then
    PIPELINE_ROOT="$(cd "${PIPELINE_ROOT}" && pwd)"
else
    PIPELINE_ROOT="$(cd "${STEP1B_BIN_DIR}/../../.." && pwd)"
fi
export PIPELINE_ROOT

STEP1B_MODULE_DIR="${PIPELINE_ROOT}/modules/step1b"
STEP1B_BIN_DIR="${STEP1B_MODULE_DIR}/bin"
STEP1B_TEMPLATE_DIR="${STEP1B_MODULE_DIR}/templates"

if [ ! -d "${STEP1B_MODULE_DIR}" ]; then
    echo "[step1b] ❌ Step1B module directory not found at ${STEP1B_MODULE_DIR}. Set PIPELINE_ROOT to the repository root." >&2
    exit 1
fi

source "${PIPELINE_ROOT}/lib/logging.sh"
source "${PIPELINE_ROOT}/lib/slurm.sh"
source "${PIPELINE_ROOT}/lib/validation.sh"
source "${STEP1B_MODULE_DIR}/lib/functions.sh"
source "${PIPELINE_ROOT}/config/pipeline_config.sh"

STEP1B_GATK_MODULE="${STEP1B_GATK_MODULE:-gatk/4.3.0.0-gcccore-11.3.0-java-11}"
LAST_STEP1B_ARRAY_MAX=0
STEP1B_FAILURE_FLAG_PATH=""
STEP1B_LAST_ERROR=""
STEP1B_RUNNING_FLAG_PATH=""

stage_step1b_reference_assets() {
    local reference_genome="$1"
    local dataset_name="$2"
    local ref_basename
    ref_basename="$(basename "${reference_genome}")"

    local shared_base="${SCRATCH_BASE_PATH%/}/${dataset_name}_shared"
    local shared_ref_dir="${shared_base}/Reference_genome"
    local shared_candidate="${shared_ref_dir}/${ref_basename}"

    mkdir -p "${shared_ref_dir}"

    if [ -f "${shared_candidate}" ]; then
        log_info "Using shared reference genome at ${shared_candidate}"
        echo "${shared_candidate}"
        return 0
    fi

    local source_dir="${PIPELINE_REFERENCE_DIR:-$(dirname "${reference_genome}")}"
    log_warn "Shared reference missing; rsync from ${source_dir}"
    if ! rsync -rhPt "${source_dir%/}/" "${shared_ref_dir}/"; then
        error_exit "Failed to copy reference genome directory from ${source_dir} to ${shared_ref_dir}. Set PIPELINE_REFERENCE_DIR correctly."
    fi

    if [ ! -f "${shared_candidate}" ]; then
        error_exit "Reference genome ${ref_basename} still missing after rsync. Verify PIPELINE_REFERENCE_FASTA."
    fi

    echo "${shared_candidate}"
}

ensure_reference_index_present() {
    local fasta_path="$1"
    local fai_path="${fasta_path}.fai"
    if [ -f "${fai_path}" ]; then
        return 0
    fi

    if ! command -v samtools >/dev/null 2>&1; then
        if [ -n "${SAMTOOLS_MODULE:-}" ]; then
            log_info "Loading samtools module ${SAMTOOLS_MODULE} for faidx generation"
            module load "${SAMTOOLS_MODULE}" >/dev/null 2>&1 || log_warn "Failed to load samtools module ${SAMTOOLS_MODULE}"
        elif command -v module >/dev/null 2>&1 && [ -n "${STEP1A_SAMTOOLS_MODULE:-}" ]; then
            log_info "Loading samtools module ${STEP1A_SAMTOOLS_MODULE} for faidx generation"
            module load "${STEP1A_SAMTOOLS_MODULE}" >/dev/null 2>&1 || log_warn "Failed to load samtools module ${STEP1A_SAMTOOLS_MODULE}"
        fi
    fi

    if command -v samtools >/dev/null 2>&1; then
        log_info "Generating FASTA index for ${fasta_path}"
        if samtools faidx "${fasta_path}"; then
            local source_dir="${PIPELINE_REFERENCE_DIR:-$(dirname "${REFERENCE_GENOME_ORIGINAL:-${fasta_path}}" )}"
            local source_fai="${source_dir%/}/$(basename "${REFERENCE_GENOME_ORIGINAL:-${fasta_path}}").fai"
            if [ ! -f "${source_fai}" ]; then
                cp -f "${fai_path}" "${source_fai}" 2>/dev/null || log_warn "Failed to copy ${fai_path} back to ${source_fai}"
            fi
            return 0
        fi
        log_warn "samtools faidx failed for ${fasta_path}; reference index still missing."
    else
        log_warn "samtools not available; cannot auto-generate ${fai_path}."
    fi

    return 1
}

step1b_failure_trap() {
    local exit_code="$1"
    if [ -z "${MASTER_LOG_FILE:-}" ] && [ -z "${STEP1B_FAILURE_FLAG_PATH:-}" ]; then
        return
    fi
    if [ "${exit_code}" -ne 0 ]; then
        local msg="${STEP1B_LAST_ERROR:-"Step 1B failed. Review logs and Data_management structure."}"
        if [ -n "${MASTER_LOG_FILE:-}" ]; then
            log_error "Step 1B failure: ${msg}"
        elif [ -n "${STEP1B_FAILURE_FLAG_PATH:-}" ]; then
            local dir
            dir="$(dirname "${STEP1B_FAILURE_FLAG_PATH}")"
            mkdir -p "${dir}"
            printf '%s\t%s\n' "$(date +%Y-%m-%dT%H:%M:%S)" "${msg}" > "${STEP1B_FAILURE_FLAG_PATH}"
        fi
        if [ -n "${STEP1B_RUNNING_FLAG_PATH:-}" ]; then
            rm -f "${STEP1B_RUNNING_FLAG_PATH}" 2>/dev/null || true
        fi
    else
        if [ -n "${STEP1B_FAILURE_FLAG_PATH:-}" ]; then
            rm -f "${STEP1B_FAILURE_FLAG_PATH}" 2>/dev/null || true
        fi
    fi
}

setup_step1b_failure_tracking() {
    local dataset="$1"
    local log_dir="${MASTER_LOG_DIR:-${LOG_BASE_PATH%/}/${dataset}}"
    mkdir -p "${log_dir}"
    STEP1B_FAILURE_FLAG_PATH="${log_dir}/step1b_failed.flag"
    rm -f "${STEP1B_FAILURE_FLAG_PATH}" 2>/dev/null || true
    export STEP1B_FAILURE_FLAG_PATH
    STEP1B_FAILURE_CONTEXT="Step1B orchestrator or array task failure"
    export STEP1B_FAILURE_CONTEXT
    STEP1B_RUNNING_FLAG_PATH="${log_dir}/step1b_running.flag"
    export STEP1B_RUNNING_FLAG_PATH
}

ensure_step1b_prereqs() {
    local resolved_gatk="${GATK_COMMAND:-}"

    if [ -n "${resolved_gatk}" ] && command -v "${resolved_gatk}" >/dev/null 2>&1; then
        log_debug "Using configured GATK command: ${resolved_gatk}"
        return 0
    fi

    if command -v gatk >/dev/null 2>&1; then
        GATK_COMMAND="$(command -v gatk)"
        export GATK_COMMAND
        log_info "Detected GATK in PATH: ${GATK_COMMAND}"
        return 0
    fi

    if ! command -v module >/dev/null 2>&1 && [ -f /etc/profile.d/modules.sh ]; then
        # shellcheck disable=SC1091
        source /etc/profile.d/modules.sh
    fi

    if command -v module >/dev/null 2>&1; then
        if module load "${STEP1B_GATK_MODULE}" >/dev/null 2>&1; then
            GATK_COMMAND="$(command -v gatk)"
            if [ -n "${GATK_COMMAND}" ]; then
                export GATK_COMMAND
                log_info "Loaded ${STEP1B_GATK_MODULE} for Step 1B: ${GATK_COMMAND}"
                return 0
            fi
        else
            log_warn "Failed to load module ${STEP1B_GATK_MODULE}"
        fi
    else
        log_warn "Environment modules command not available; cannot auto-load ${STEP1B_GATK_MODULE}"
    fi

    error_exit "GATK command not found. Load ${STEP1B_GATK_MODULE} manually or set GATK_COMMAND."
}

main() {
    local dataset_name="$1"
    local rdm_base_path="$2"
    STEP1B_LAST_ERROR="Step 1B failed before reason was recorded."
    setup_step1b_failure_tracking "${dataset_name}"
    trap 'step1b_failure_trap $?' EXIT

    if [ -z "${MASTER_LOG_DIR:-}" ]; then
    init_logging "step1b" "pipeline"
    else
        LOG_TO_FILE="false"
        LOG_TO_CONSOLE="true"
        MODULE_NAME="step1b"
        LOG_LEVEL="${LOG_LEVEL:-INFO}"
        log_info "Routing Step 1B logs to master log."
    fi
    ensure_step1b_prereqs

    log_info "Starting Step 1B execution"
    log_info "Dataset: ${dataset_name}"
    log_info "RDM Base Path: ${rdm_base_path}"

    validate_step1b_inputs "${rdm_base_path}"

    local reference_genome
    reference_genome="$(get_reference_fasta)"
    if [ ! -f "${reference_genome}" ]; then
        log_error "Reference genome not found at ${reference_genome}"
        STEP1B_LAST_ERROR="Reference genome not found: ${reference_genome}"
        exit 1
    fi

    local staged_reference
    staged_reference="$(stage_step1b_reference_assets "${reference_genome}" "${dataset_name}")"
    export REFERENCE_GENOME="${staged_reference}"
    export PIPELINE_REFERENCE_FASTA="${staged_reference}"
    export REFERENCE_GENOME_INDEX="${staged_reference}.fai"
    export REFERENCE_GENOME_DICT="${staged_reference%.*}.dict"
    if ! ensure_reference_index_present "${staged_reference}"; then
        STEP1B_LAST_ERROR="Reference genome index missing for ${staged_reference} (expected ${staged_reference}.fai)."
        exit 1
    fi
    reference_genome="${staged_reference}"

    local chromosome_list
    if [ -n "${STEP1B_CHROMOSOME_LIST:-}" ]; then
        log_info "Using STEP1B_CHROMOSOME_LIST override: ${STEP1B_CHROMOSOME_LIST}"
        chromosome_list="$(cat "${STEP1B_CHROMOSOME_LIST}")"
    else
    chromosome_list="$(get_chromosome_list "${reference_genome}")"
    fi
    local chromosome_count
    chromosome_count=$(printf '%s\n' "${chromosome_list}" | grep -c '.')
    log_info "Reference genome: ${reference_genome}"
    log_info "Chromosome list lines: ${chromosome_count}"
    if [ "${chromosome_count}" -le 0 ]; then
        log_error "No chromosomes detected for Step 1B. Ensure ${reference_genome} exists and has an index (samtools faidx)."
        STEP1B_LAST_ERROR="Chromosome discovery returned 0 entries for ${reference_genome}. Confirm reference/index per Data_management inventory."
        exit 1
    fi

    log_info "Chromosomes to process: ${chromosome_count}"

    local step1b_status
    step1b_status=$(check_step1b_status "${rdm_base_path}")
    case "${step1b_status}" in
        "Complete")
            log_info "Step 1B outputs already present."
            if ! confirm_action "Outputs detected. Do you want to rerun Step 1B?"; then
                log_info "Step 1B execution cancelled."
                exit 0
            fi
            ;;
        "Partial"|"Incomplete"|"Not Started")
            log_info "Step 1B will run (status: ${step1b_status})"
            ;;
    esac

    local scratch_chr_dir="${PIPELINE_WORK_DIR%/}/step1b/${dataset_name}"
    mkdir -p "${scratch_chr_dir}"

    local chromosome_file="${scratch_chr_dir}/chromosomes_${dataset_name}_$(date +%Y%m%d_%H%M%S).txt"
    printf '%s\n' "${chromosome_list}" > "${chromosome_file}"

    local slurm_script
    if ! slurm_script=$(create_step1b_slurm_script "${dataset_name}" "${chromosome_count}"); then
        log_error "Failed to generate Step 1B SLURM script"
        STEP1B_LAST_ERROR="Failed to generate Step 1B SLURM script."
        exit 1
    fi
    if [ ! -f "${slurm_script}" ]; then
        log_error "Generated Step 1B SLURM script not found: ${slurm_script}"
        STEP1B_LAST_ERROR="Generated Step 1B SLURM script missing: ${slurm_script}"
        exit 1
    fi

    local job_id
    job_id=$(submit_job "${slurm_script}" "${rdm_base_path} ${chromosome_file}" "${dataset_name}" "1B")

    if [ $? -eq 0 ]; then
        log_slurm_submission "${job_id}" "${slurm_script}" "${dataset_name}" "${chromosome_count}"
        local array_range_max="${LAST_STEP1B_ARRAY_MAX:-$((chromosome_count - 1))}"

        if [ -n "${STEP1B_RUNNING_FLAG_PATH:-}" ]; then
            mkdir -p "$(dirname "${STEP1B_RUNNING_FLAG_PATH}")"
            printf '%s\tJobID=%s\tChromosomes=%s\n' "$(date +%Y-%m-%dT%H:%M:%S)" "${job_id}" "${chromosome_count}" > "${STEP1B_RUNNING_FLAG_PATH}" 2>/dev/null || true
        fi

        echo "✓ Step 1B submitted successfully!"
        echo "Job ID: ${job_id}"
        echo "Array jobs: 0-${array_range_max}"
        echo "Monitor with: squeue -u \$USER"
        echo "Cancel with: scancel ${job_id}"
        local log_hint="${MASTER_LOG_DIR:-${LOG_BASE_PATH%/}/${dataset_name}}"
        echo "View logs: ls ${log_hint}"
    else
        log_error "Failed to submit Step 1B."
        STEP1B_LAST_ERROR="Failed to submit Step 1B orchestrator via sbatch."
        exit 1
    fi
}

execute_step1b_pipeline() {
    local chromosome="$1"
    local rdm_base_path="$2"
    local dataset_name="${3:-$(basename "${rdm_base_path}")}"

    log_info "Processing chromosome ${chromosome}"

    local reference_genome
    reference_genome="$(get_reference_fasta)"
    if [ ! -f "${reference_genome}" ]; then
        error_exit "Reference genome not found at ${reference_genome}"
    fi

    local dataset_folder
    dataset_folder="$(basename "${rdm_base_path}")"
    local workdir="${SCRATCH_BASE_PATH%/}/${dataset_folder}_1B"

    ensure_step1b_workdir "${workdir}"

    # Copy reference genome to working directory (matching original script approach)
    # Uses symlinks if shared reference genome exists (for array jobs)
    local local_reference_genome
    local_reference_genome="$(copy_reference_genome_to_workdir "${reference_genome}" "${workdir}" "${dataset_name}")"

    local sample_map="${workdir}/sample_map.txt"
    local lock_dir="${sample_map}.lock"

    if [ ! -s "${sample_map}" ]; then
        if mkdir "${lock_dir}" 2>/dev/null; then
            build_sample_map "${sample_map}" "${rdm_base_path}"
            rmdir "${lock_dir}"
        else
            log_info "Waiting for sample map to be generated..."
            while [ ! -s "${sample_map}" ]; do
                sleep 5
            done
        fi
    fi

    local memory="${STEP1B_MEMORY:-36G}"

    run_genomics_db_import "${chromosome}" "${workdir}" "${local_reference_genome}" "${sample_map}" "${memory}"

    local output_file="${workdir}/${chromosome}_consolidated.vcf.gz"
    run_genotype_gvcfs "${chromosome}" "${workdir}" "${local_reference_genome}" "${memory}" "${output_file}"

    copy_consolidated_vcf "${output_file}" "${rdm_base_path}" "${chromosome}"

    cleanup_chromosome_workspace "${workdir}" "${chromosome}"

    log_info "Chromosome ${chromosome} processing complete."
}

validate_step1b_inputs() {
    local rdm_base_path="$1"

    log_info "Validating Step 1B inputs"

    if [ ! -d "${rdm_base_path}" ]; then
        log_error "RDM base path does not exist: ${rdm_base_path}"
        STEP1B_LAST_ERROR="RDM base path missing: ${rdm_base_path}"
        exit 1
    fi

    local gvcf_dir="${rdm_base_path}/5.Individual_VCF"
    if [ ! -d "${gvcf_dir}" ]; then
        log_error "Directory not found: ${gvcf_dir}"
        STEP1B_LAST_ERROR="Missing 5.Individual_VCF under ${rdm_base_path}."
        exit 1
    fi

    local gvcf_count
    gvcf_count=$(find "${gvcf_dir}" -name "*_raw.g.vcf.gz" -type f | wc -l)
    if [ "${gvcf_count}" -eq 0 ]; then
        log_error "No *_raw.g.vcf.gz files found in ${gvcf_dir}"
        STEP1B_LAST_ERROR="No raw GVCFs found in ${gvcf_dir}."
        exit 1
    fi

    log_info "Detected ${gvcf_count} GVCF files."
}

create_step1b_slurm_script() {
    local dataset_name="$1"
    local chromosome_count="$2"

    local config
    config=$(get_step1b_config)
    declare -A config_map=()
    while IFS='=' read -r key value; do
        [ -z "${key}" ] && continue
        config_map["${key}"]="${value}"
    done <<< "${config}"

    local array_max=$((chromosome_count - 1))
    local configured_array_limit="${config_map[ARRAY_MAX]:-}"
    if [[ -n "${configured_array_limit}" && "${configured_array_limit}" -gt 0 && "${array_max}" -ge "${configured_array_limit}" ]]; then
        log_warn "Chromosome count (${chromosome_count}) exceeds configured array limit (${configured_array_limit}); truncating array."
        array_max=$((configured_array_limit - 1))
    fi
    LAST_STEP1B_ARRAY_MAX="${array_max}"

    mkdir -p "${PIPELINE_SLURM_SCRIPT_DIR}"
    local slurm_script="${PIPELINE_SLURM_SCRIPT_DIR}/Apple_GATK_1B_${dataset_name}_$(date +%Y%m%d_%H%M%S).sh"
    local template="${STEP1B_TEMPLATE_DIR}/step1b_array.sh"

    local config_string="job_name=Apple_GATK_1B_${dataset_name}
account=${config_map[ACCOUNT]}
partition=${config_map[PARTITION]}
nodes=${config_map[NODES]}
ntasks=${config_map[NTASKS]}
cpus_per_task=${config_map[CPUS]}
time_limit=${config_map[TIME]}
memory=${config_map[MEMORY]}
array_max=${array_max}"
    if [ -n "${config_map[QOS]:-}" ]; then
        config_string="${config_string}
qos=${config_map[QOS]}"
    fi

    create_slurm_script "${template}" "${config_string}" "${slurm_script}" "${dataset_name}" "1B" "${PIPELINE_ROOT}"
    echo "${slurm_script}"
}

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    main "$@"
fi
