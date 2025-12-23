#!/bin/bash
# =============================================================================
# STEP 1C - BEAGLE IMPUTATION
# =============================================================================

STEP1C_BIN_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STEP1C_MODULE_DIR="$(cd "${STEP1C_BIN_DIR}/.." && pwd)"
STEP1C_TEMPLATE_DIR="${STEP1C_MODULE_DIR}/templates"
if [ -n "${PIPELINE_ROOT:-}" ]; then
    if [ -d "${PIPELINE_ROOT}" ]; then
        PIPELINE_ROOT="$(cd "${PIPELINE_ROOT}" && pwd)"
    else
        echo "[step1c] ⚠️  Provided PIPELINE_ROOT (${PIPELINE_ROOT}) does not exist; falling back to module-relative path." >&2
        PIPELINE_ROOT="$(cd "${STEP1C_BIN_DIR}/../../.." && pwd)"
    fi
else
    PIPELINE_ROOT="$(cd "${STEP1C_BIN_DIR}/../../.." && pwd)"
fi
export PIPELINE_ROOT

source "${PIPELINE_ROOT}/lib/logging.sh"
source "${PIPELINE_ROOT}/lib/slurm.sh"
source "${STEP1C_MODULE_DIR}/lib/functions.sh"
source "${PIPELINE_ROOT}/config/pipeline_config.sh"

LAST_STEP1C_ARRAY_MAX=0

main() {
    local dataset_name="$1"
    local rdm_base_path="$2"
    local gene_map_file="${GENE_MAP_FILE:-}"
    local output_dir_override="${STEP1C_OUTPUT_DIR:-}"

    init_logging "step1c" "pipeline" "${dataset_name}"

    if [ -z "${dataset_name}" ] || [ -z "${rdm_base_path}" ]; then
        log_error "Usage: step1c main <dataset_name> <rdm_base_path>"
        exit 1
    fi

    validate_step1c_inputs "${rdm_base_path}"

    local reference_genome
    reference_genome="$(get_reference_fasta)"
    if [ ! -f "${reference_genome}" ]; then
        log_error "Reference genome not found: ${reference_genome}"
        exit 1
    fi

    local output_dir
    if [ -n "${output_dir_override}" ]; then
        output_dir="${output_dir_override}"
    else
        output_dir="$(default_step1c_output_dir "${rdm_base_path}")"
    fi

    mkdir -p "${output_dir}"

    local scratch_manifest_dir="${PIPELINE_WORK_DIR%/}/step1c/${dataset_name}"
    mkdir -p "${scratch_manifest_dir}"

    local manifest_file="${scratch_manifest_dir}/vcf_manifest_${dataset_name}_$(date +%Y%m%d_%H%M%S).txt"
    create_vcf_manifest "${rdm_base_path}" "${manifest_file}"
    local manifest_count
    manifest_count=$(grep -c . "${manifest_file}")
    log_info "VCF manifest entries: ${manifest_count}"

    local slurm_script
    if ! slurm_script=$(create_step1c_slurm_script "${dataset_name}" "${manifest_count}"); then
        log_error "Failed to generate Step 1C SLURM script"
        exit 1
    fi
    if [ ! -f "${slurm_script}" ]; then
        log_error "Generated Step 1C SLURM script not found: ${slurm_script}"
        exit 1
    fi

    local thread_count="${STEP1C_CPUS_PER_TASK:-8}"
    local memory_string="${STEP1C_MEMORY:-48G}"
    local memory_gb="${memory_string%G}"
    local self_impute="${STEP1C_SELF_IMPUTE:-false}"
    local self_impute_lc
    self_impute_lc="$(echo "${self_impute}" | tr '[:upper:]' '[:lower:]')"
    local impute_flag
    case "${self_impute_lc}" in
        true|1|yes|y)
            impute_flag="true"
            ;;
        false|0|no|n|"")
            impute_flag="false"
            ;;
        *)
            log_warn "Unrecognized STEP1C_SELF_IMPUTE=${self_impute}; defaulting to false (phasing-only)."
            impute_flag="false"
            ;;
    esac
    log_info "Step 1C Beagle mode: self-impute=${impute_flag}"

    # Preserve positional args even when gene map is omitted
    local gene_map_arg="${gene_map_file:-__NO_GENE_MAP__}"

    local cmd_args=("${rdm_base_path}" "${manifest_file}" "${output_dir}" "${reference_genome}"
                    "${gene_map_arg}" "${dataset_name}_beagle" "${thread_count}"
                    "${memory_gb}" "${impute_flag}")

    local job_id
    job_id=$(submit_job "${slurm_script}" "${cmd_args[*]}" "${dataset_name}" "1C")

    if [ $? -eq 0 ]; then
        log_info "Step 1C job submitted successfully (ID: ${job_id})"
        log_info "Array range: 0-${LAST_STEP1C_ARRAY_MAX}"
    else
        log_error "Failed to submit Step 1C job"
        exit 1
    fi
}

create_step1c_slurm_script() {
    local dataset_name="$1"
    local manifest_count="$2"

    if [ -z "${manifest_count}" ] || [ "${manifest_count}" -le 0 ]; then
        log_error "Manifest is empty; cannot create Step 1C array job"
        return 1
    fi

    local config
    config=$(get_step1c_config)
    declare -A config_map=()
    while IFS='=' read -r key value; do
        [ -z "${key}" ] && continue
        config_map["${key}"]="${value}"
    done <<< "${config}"

    local memory="${config_map[MEMORY]}"
    local cpus="${config_map[CPUS]}"
    local time_limit="${config_map[TIME]}"
    local array_max=$((manifest_count - 1))
    local configured_array_limit="${config_map[ARRAY_MAX]:-}"
    if [[ -n "${configured_array_limit}" && "${configured_array_limit}" -gt 0 && "${array_max}" -ge "${configured_array_limit}" ]]; then
        log_warn "VCF count (${manifest_count}) exceeds configured array limit (${configured_array_limit}); truncating array."
        array_max=$((configured_array_limit - 1))
    fi
    LAST_STEP1C_ARRAY_MAX="${array_max}"

    mkdir -p "${PIPELINE_SLURM_SCRIPT_DIR}"
    local slurm_script="${PIPELINE_SLURM_SCRIPT_DIR}/Apple_GATK_1C_${dataset_name}_$(date +%Y%m%d_%H%M%S).sh"
    local template="${STEP1C_TEMPLATE_DIR}/step1c_job.sh"

    if [ ! -f "${template}" ]; then
        log_error "SLURM template not found: ${template}"
        return 1
    fi
    log_info "Using Step 1C template: ${template}"

    local config_string="job_name=Apple_GATK_1C_${dataset_name}
account=${config_map[ACCOUNT]}
partition=${config_map[PARTITION]}
nodes=${config_map[NODES]}
ntasks=${config_map[NTASKS]}
cpus_per_task=${cpus}
time_limit=${time_limit}
memory=${memory}
array_max=${array_max}"
    if [ -n "${config_map[QOS]:-}" ]; then
        config_string="${config_string}
qos=${config_map[QOS]}"
    fi
    if [ -n "${config_map[CONSTRAINT]:-}" ]; then
        config_string="${config_string}
constraint=${config_map[CONSTRAINT]}"
    fi

    create_slurm_script "${template}" "${config_string}" "${slurm_script}" "${dataset_name}" "1C" "${PIPELINE_ROOT}"
    echo "${slurm_script}"
}

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    main "$@"
fi
