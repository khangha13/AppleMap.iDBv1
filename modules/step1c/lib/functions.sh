#!/bin/bash
# =============================================================================
# STEP 1C FUNCTIONS - BEAGLE IMPUTATION
# =============================================================================

# Validate Step 1C inputs
validate_step1c_inputs() {
    local rdm_base_path="$1"

    if [ ! -d "${rdm_base_path}" ]; then
        log_error "RDM base path does not exist: ${rdm_base_path}"
        exit 1
    fi

    local joint_dir="${rdm_base_path}/7.Consolidated_VCF"
    if [ ! -d "${joint_dir}" ]; then
        log_error "Joint VCF directory not found: ${joint_dir}"
        exit 1
    fi

    local joint_count
    joint_count=$(find "${joint_dir}" -maxdepth 1 -name '*.vcf.gz' -type f | wc -l)
    if [ "${joint_count}" -eq 0 ]; then
        log_error "No *.vcf.gz files found in ${joint_dir}"
        exit 1
    fi

    log_info "Detected ${joint_count} joint VCF files for imputation."
}

# Create a file containing absolute paths to VCFs for Beagle
create_vcf_manifest() {
    local rdm_base_path="$1"
    local manifest_path="$2"

    local joint_dir="${rdm_base_path}/7.Consolidated_VCF"
    # Include only Chr01–Chr17; allow suffixes like _consolidated before .vcf.gz
    find "${joint_dir}" -maxdepth 1 -type f -name 'Chr*.vcf.gz' \
        | awk -F/ '{
            base=$NF;
            if (base ~ /^Chr(0[1-9]|1[0-7]).*\.vcf\.gz$/) print;
        }' \
        | sort > "${manifest_path}"

    if [ ! -s "${manifest_path}" ]; then
        error_exit "VCF manifest ${manifest_path} is empty (no Chr01–Chr17 VCFs found)."
    fi
}

# Determine default output directory
default_step1c_output_dir() {
    local rdm_base_path="$1"
    echo "${rdm_base_path}/8.Imputated_VCF_BEAGLE"
}
