#!/bin/bash
# =============================================================================
# STEP 1D SUBMISSION WRAPPER
# =============================================================================

STEP1D_BIN_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODULE_DIR="$(cd "${STEP1D_BIN_DIR}/.." && pwd)"
if [ -n "${PIPELINE_ROOT:-}" ]; then
    if [ -d "${PIPELINE_ROOT}" ]; then
        PIPELINE_ROOT="$(cd "${PIPELINE_ROOT}" && pwd)"
    else
        echo "[step1d] ⚠️  Provided PIPELINE_ROOT (${PIPELINE_ROOT}) does not exist; falling back to module-relative path." >&2
        PIPELINE_ROOT="$(cd "${STEP1D_BIN_DIR}/../../.." && pwd)"
    fi
else
    PIPELINE_ROOT="$(cd "${STEP1D_BIN_DIR}/../../.." && pwd)"
fi
export PIPELINE_ROOT
SCRIPT_DIR="${STEP1D_BIN_DIR}"

source "${PIPELINE_ROOT}/lib/logging.sh"
source "${PIPELINE_ROOT}/lib/slurm.sh"
source "${PIPELINE_ROOT}/config/pipeline_config.sh"

main() {
    local dataset_name="$1"
    local vcf_dir="$2"
    shift 2 || true

    local beagle_flag=false
    local dry_run_flag=false
    local remove_rel_flag=false
    local mode="qc"
    local mode_set=false
    while (( "$#" )); do
        case "$1" in
            --beagle)
                beagle_flag=true
                ;;
            --dry-run|-n)
                dry_run_flag=true
                ;;
            --qc)
                if ${mode_set}; then
                    log_error "Multiple modes provided; choose one of --qc, --PCA, or --duplicate-check."
                    exit 1
                fi
                mode="qc"
                mode_set=true
                ;;
            --PCA|--pca)
                if ${mode_set}; then
                    log_error "Multiple modes provided; choose one of --qc, --PCA, or --duplicate-check."
                    exit 1
                fi
                mode="pca"
                mode_set=true
                ;;
            --duplicate-check)
                if ${mode_set}; then
                    log_error "Multiple modes provided; choose one of --qc, --PCA, or --duplicate-check."
                    exit 1
                fi
                mode="duplicate-check"
                mode_set=true
                ;;
            --remove-relatives)
                remove_rel_flag=true
                ;;
            *)
                log_warn "Unknown option ignored: $1"
                ;;
        esac
        shift || true
    done

    if ${remove_rel_flag} && [ "${mode}" != "pca" ]; then
        log_error "--remove-relatives requires --PCA."
        exit 1
    fi

    init_logging "step1d" "pipeline" "${dataset_name}"

    if [ -z "${dataset_name}" ] || [ -z "${vcf_dir}" ]; then
        log_error "Usage: step1d main <dataset_name> <vcf_directory> [--beagle] [--dry-run] [--qc|--PCA|--duplicate-check] [--remove-relatives]"
        exit 1
    fi

    if [ ! -d "${vcf_dir}" ]; then
        log_error "VCF directory not found: ${vcf_dir}"
        exit 1
    fi

    export VCF_DIR="${vcf_dir}"
    export WORK_DIR="${WORK_DIR_OVERRIDE:-${vcf_dir}}"
    export R_SCRIPTS_DIR="${R_SCRIPTS_DIR_OVERRIDE:-${MODULE_DIR}/Rscripts}"

    local config
    config=$(get_step1d_config)
    declare -A config_map=()
    while IFS='=' read -r key value; do
        [ -z "${key}" ] && continue
        config_map["${key}"]="${value}"
    done <<< "${config}"

    local slurm_script
    if ! slurm_script=$(create_step1d_slurm_script "${dataset_name}"); then
        log_error "Failed to generate Step 1D SLURM script"
        exit 1
    fi
    if [ ! -f "${slurm_script}" ]; then
        log_error "Generated Step 1D SLURM script not found: ${slurm_script}"
        exit 1
    fi

    local job_args=()
    $beagle_flag && job_args+=("--beagle")
    $dry_run_flag && job_args+=("--dry-run")
    job_args+=("--${mode}")
    $remove_rel_flag && job_args+=("--remove-relatives")

    local job_id
    job_id=$(submit_job "${slurm_script}" "${job_args[*]}" "${dataset_name}" "1D")

    if [ $? -eq 0 ]; then
        log_info "Step 1D job submitted successfully (ID: ${job_id})"
    else
        log_error "Failed to submit Step 1D job"
        exit 1
    fi
}

create_step1d_slurm_script() {
    local dataset_name="$1"

    local config
    config=$(get_step1d_config)
    declare -A config_map=()
    while IFS='=' read -r key value; do
        [ -z "${key}" ] && continue
        config_map["${key}"]="${value}"
    done <<< "${config}"

    mkdir -p "${PIPELINE_SLURM_SCRIPT_DIR}"
    local slurm_script="${PIPELINE_SLURM_SCRIPT_DIR}/Apple_GATK_1D_${dataset_name}_$(date +%Y%m%d_%H%M%S).sh"
    local template="${SCRIPT_DIR}/../templates/master_vcf_analysis.sh"

    local config_string="job_name=Apple_GATK_1D_${dataset_name}
account=${config_map[ACCOUNT]}
partition=${config_map[PARTITION]}
nodes=${config_map[NODES]}
ntasks=${config_map[NTASKS]}
cpus_per_task=${config_map[CPUS]}
time_limit=${config_map[TIME]}
memory=${config_map[MEMORY]}
array_max=0"
    if [ -n "${config_map[QOS]:-}" ]; then
        config_string="${config_string}
qos=${config_map[QOS]}"
    fi
    if [ -n "${config_map[CONSTRAINT]:-}" ]; then
        config_string="${config_string}
constraint=${config_map[CONSTRAINT]}"
    fi

    create_slurm_script "${template}" "${config_string}" "${slurm_script}" "${dataset_name}" "1D" "${PIPELINE_ROOT}"
    echo "${slurm_script}"
}

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    main "$@"
fi
