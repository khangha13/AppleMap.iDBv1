#!/bin/bash
# =============================================================================
# VALIDATION FUNCTIONS FOR GATK PIPELINE
# =============================================================================
# Common validation functions used across all pipeline steps

# Function to validate directory structure
validate_directory_structure() {
    local base_path=$1
    local required_dirs=("1.FASTQ" "2.FASTQC_pre_trimmed" "3.FASTQC_post_trimmed" "4.BAM" "5.Individual_VCF" "6.genomicsdb_output" "7.Consolidated_VCF" "8.Imputated_VCF_BEAGLE" "9.Imputation_QC")
    
    print_message $BLUE "Validating directory structure for: $base_path"
    
    for dir in "${required_dirs[@]}"; do
        if [ ! -d "$base_path/$dir" ]; then
            print_message $RED "✗ Missing directory: $dir"
        else
            print_message $GREEN "✓ Directory exists: $dir"
        fi
    done
    
    # Check if 1.FASTQ contains sample files
    if [ -d "$base_path/1.FASTQ" ]; then
        local sample_count=$(find "$base_path/1.FASTQ" -name "*_1.fastq.gz" | wc -l)
        if [ $sample_count -eq 0 ]; then
            print_message $RED "Warning: No paired-end FASTQ files found in $base_path/1.FASTQ"
            print_message $YELLOW "Expected format: SAMPLE_1.fastq.gz and SAMPLE_2.fastq.gz"
        else
            print_message $GREEN "✓ Found $sample_count paired-end sample files"
        fi
    fi
}

# Function to validate sample files
validate_sample_files() {
    local fastq_dir=$1
    
    print_message $BLUE "Scanning for sample files in: $fastq_dir"
    
    # Find all _1.fastq.gz files and extract sample names
    local samples=($(find "$fastq_dir" -name "*_1.fastq.gz" -exec basename {} \; | sed 's/_1\.fastq\.gz$//' | sort))
    
    if [ ${#samples[@]} -eq 0 ]; then
        print_message $RED "Error: No sample files found in $fastq_dir"
        print_message $YELLOW "Expected format: SAMPLE_1.fastq.gz and SAMPLE_2.fastq.gz"
        return 1
    fi
    
    # Verify paired-end files exist for each sample
    local valid_samples=()
    for sample in "${samples[@]}"; do
        if [ -f "$fastq_dir/${sample}_1.fastq.gz" ] && [ -f "$fastq_dir/${sample}_2.fastq.gz" ]; then
            valid_samples+=("$sample")
            print_message $GREEN "✓ Valid sample: $sample"
        else
            print_message $RED "✗ Incomplete sample: $sample (missing paired-end file)"
        fi
    done
    
    if [ ${#valid_samples[@]} -eq 0 ]; then
        print_message $RED "Error: No valid paired-end samples found"
        return 1
    fi
    
    print_message $GREEN "✓ Found ${#valid_samples[@]} valid samples"
    return 0
}

# Function to check Step 1A completion status
check_step1a_status() {
    local rdm_base_path=$1
    
    log_info "Checking Step 1A completion status..."
    
    # Normalize path (trim CRs and trailing whitespace) and remove trailing slash
    rdm_base_path="$(printf '%s' "$rdm_base_path" | tr -d '\r' | sed -e 's/[[:space:]]\+$//')"
    local individual_vcf_dir="${rdm_base_path%/}/5.Individual_VCF"
    
    if [ ! -d "$individual_vcf_dir" ]; then
        log_warn "Step 1A directory does not exist: $individual_vcf_dir"
        echo "Not Started"
        return
    fi
    
    # Check for individual VCF files
    local vcf_files=($(find "$individual_vcf_dir" -name "*_genotyped.vcf.gz" -type f 2>/dev/null))
    local gvcf_files=($(find "$individual_vcf_dir" -name "*_raw.g.vcf.gz" -type f 2>/dev/null))
    
    local vcf_count=${#vcf_files[@]}
    local gvcf_count=${#gvcf_files[@]}
    
    if [ $vcf_count -gt 0 ]; then
        # Check if all VCF files have index files
        local complete_vcfs=0
        for vcf_file in "${vcf_files[@]}"; do
            if [ -f "${vcf_file}.tbi" ]; then
                complete_vcfs=$((complete_vcfs + 1))
            fi
        done
        
        if [ $complete_vcfs -eq $vcf_count ]; then
            log_info "Step 1A complete: $complete_vcfs genotyped VCF files"
            echo "Complete"
        else
            log_warn "Step 1A partially complete: $complete_vcfs/$vcf_count VCF files"
            echo "Partial"
        fi
    elif [ $gvcf_count -gt 0 ]; then
        log_warn "Step 1A partially complete: $gvcf_count GVCF files (missing genotyped VCFs)"
        echo "Partial"
    else
        log_info "Step 1A not started: no VCF files found"
        echo "Not Started"
    fi
}

# Return list of incomplete samples (one per line) given an RDM base path
get_incomplete_samples() {
    local rdm_base_path="$1"
    rdm_base_path="$(printf '%s' "$rdm_base_path" | tr -d '\r' | sed -e 's/[[:space:]]\+$//')"
    local fastq_dir="${rdm_base_path%/}/1.FASTQ"
    local indiv_dir="${rdm_base_path%/}/5.Individual_VCF"
    local tmp_list
    tmp_list="$(mktemp 2>/dev/null || echo "/tmp/sample_list.$$")"
    local -a samples=()

    if [ -d "$fastq_dir" ]; then
        # Build sample list from FASTQ directory
        if create_sample_list "$fastq_dir" "$tmp_list" >/dev/null 2>&1; then
            mapfile -t samples < "$tmp_list"
        fi
        rm -f "$tmp_list" 2>/dev/null || true
    fi

    # Fallback: infer from VCF/GVCF basenames if no samples yet
    if [ ${#samples[@]} -eq 0 ] && [ -d "$indiv_dir" ]; then
        # Collect unique basenames from either raw.g.vcf.gz or genotyped.vcf.gz
        mapfile -t samples < <(find "$indiv_dir" -maxdepth 1 -type f \( -name "*_raw.g.vcf.gz" -o -name "*_genotyped.vcf.gz" \) -printf "%f\n" 2>/dev/null | sed -E 's/(_raw\.g|_genotyped)\.vcf\.gz$//' | sort -u)
    fi

    # If still none, nothing to report
    if [ ${#samples[@]} -eq 0 ]; then
        return 0
    fi

    local sample
    for sample in "${samples[@]}"; do
        local vcf="${indiv_dir}/${sample}_genotyped.vcf.gz"
        if [ -f "$vcf" ] && [ -f "${vcf}.tbi" ]; then
            continue
        fi
        echo "$sample"
    done
}

# Function to check Step 1B completion status
check_step1b_status() {
    local rdm_base_path=$1
    
    log_info "Checking Step 1B completion status..."
    
    local dataset_name
    dataset_name="$(basename "${rdm_base_path%/}")"
    local failure_flag="${LOG_BASE_PATH%/}/${dataset_name}/step1b_failed.flag"
    if [ -f "${failure_flag}" ]; then
        local reason
        reason="$(cat "${failure_flag}" 2>/dev/null || true)"
        log_error "Step 1B failure flag detected${reason:+: ${reason}}"
        echo "Failed"
        return
    fi

    local consolidated_dir="${rdm_base_path}/7.Consolidated_VCF"
    
    if [ ! -d "$consolidated_dir" ]; then
        log_warn "Step 1B directory does not exist: $consolidated_dir"
        echo "Not Started"
        return
    fi
    
    # Check for consolidated VCF files (chromosome-based)
    local vcf_files=($(find "$consolidated_dir" -name "*.vcf.gz" -type f))
    local vcf_count=${#vcf_files[@]}
    
    if [ $vcf_count -eq 0 ]; then
        log_warn "No consolidated VCF files found in Step 1B directory"
        echo "Not Started"
        return
    fi
    
    # Check if files have corresponding index files
    local complete_files=0
    for vcf_file in "${vcf_files[@]}"; do
        if [ -f "${vcf_file}.tbi" ]; then
            complete_files=$((complete_files + 1))
        fi
    done
    
    if [ $complete_files -eq $vcf_count ]; then
        log_info "Step 1B complete: $complete_files consolidated VCF files"
        echo "Complete"
    elif [ $complete_files -gt 0 ]; then
        log_warn "Step 1B partially complete: $complete_files/$vcf_count VCF files"
        echo "Partial"
    else
        log_warn "Step 1B incomplete: no complete VCF files found"
        echo "Incomplete"
    fi
}

# Function to check Step 1C completion status
check_step1c_status() {
    local rdm_base_path=$1

    log_info "Checking Step 1C (Beagle) completion status..."

    local imputed_dir="${rdm_base_path}/8.Imputated_VCF_BEAGLE"

    if [ ! -d "${imputed_dir}" ]; then
        log_warn "Step 1C output directory does not exist: ${imputed_dir}"
        echo "Not Started"
        return
    fi

    local vcf_files=($(find "${imputed_dir}" -maxdepth 1 -name '*.vcf.gz' -type f))
    local vcf_count=${#vcf_files[@]}

    if [ "${vcf_count}" -eq 0 ]; then
        log_warn "No Beagle-imputed VCF files found in ${imputed_dir}"
        echo "Not Started"
        return
    fi

    local complete_files=0
    for vcf_file in "${vcf_files[@]}"; do
        if [ -f "${vcf_file}.tbi" ]; then
            complete_files=$((complete_files + 1))
        fi
    done

    if [ "${complete_files}" -eq "${vcf_count}" ]; then
        log_info "Step 1C complete: ${complete_files} imputed VCF files"
        echo "Complete"
    elif [ "${complete_files}" -gt 0 ]; then
        log_warn "Step 1C partially complete: ${complete_files}/${vcf_count} VCFs indexed"
        echo "Partial"
    else
        log_warn "Step 1C incomplete: indexes missing for all VCFs"
        echo "Incomplete"
    fi
}

# Function to check Step 1D completion status
check_step1d_status() {
    local rdm_base_path=$1

    log_info "Checking Step 1D (QC analysis) completion status..."

    local metrics_dir="${rdm_base_path}/9.Imputation_QC"

    if [ ! -d "${metrics_dir}" ]; then
        log_warn "Step 1D metrics directory does not exist: ${metrics_dir}"
        echo "Not Started"
        return
    fi

    local metrics_file="${metrics_dir}/variant_site_metrics.tsv"
    local depth_plots_dir="${metrics_dir}/depth_plots"
    local missing_plots_dir="${metrics_dir}/missingness_plots"

    local has_metrics=0
    local has_plots=0

    if [ -f "${metrics_file}" ]; then
        has_metrics=1
    fi

    if [ -d "${depth_plots_dir}" ] && [ -d "${missing_plots_dir}" ]; then
        has_plots=1
    fi

    if [ "${has_metrics}" -eq 1 ] && [ "${has_plots}" -eq 1 ]; then
        log_info "Step 1D complete: metrics and plots detected in ${metrics_dir}"
        echo "Complete"
    elif [ "${has_metrics}" -eq 1 ] || [ "${has_plots}" -eq 1 ]; then
        log_warn "Step 1D partially complete: some outputs detected"
        echo "Partial"
    else
        log_warn "Step 1D not started: no outputs detected"
        echo "Not Started"
    fi
}

# Function to check pipeline completion status for samples
check_pipeline_completion() {
    local rdm_base_path=$1
    local sample_list_file=$2
    
    log_info "Checking pipeline completion status..."
    
    local completed_samples=()
    local incomplete_samples=()
    local total_samples=0
    
    # Read sample list
    while IFS= read -r sample; do
        [ -z "$sample" ] && continue  # Skip empty lines
        total_samples=$((total_samples + 1))
        
        # Check for final output files (genotyped VCF)
        local genotyped_vcf="${rdm_base_path}/5.Individual_VCF/${sample}_genotyped.vcf.gz"
        local raw_gvcf="${rdm_base_path}/5.Individual_VCF/${sample}_raw.g.vcf.gz"
        
        if [ -f "$genotyped_vcf" ] && [ -f "${genotyped_vcf}.tbi" ]; then
            completed_samples+=("$sample")
            log_debug "Sample completed: $sample"
        elif [ -f "$raw_gvcf" ] && [ -f "${raw_gvcf}.tbi" ]; then
            # Has GVCF but not genotyped VCF - partially complete
            incomplete_samples+=("$sample")
            log_warn "Sample partially complete: $sample (has GVCF, missing genotyped VCF)"
        else
            incomplete_samples+=("$sample")
            log_debug "Sample incomplete: $sample"
        fi
    done < "$sample_list_file"
    
    log_info "Pipeline completion summary - Completed: ${#completed_samples[@]}, Incomplete: ${#incomplete_samples[@]}, Total: $total_samples"
    
    # Return completion status
    if [ ${#completed_samples[@]} -eq $total_samples ]; then
        echo "all_complete"
    elif [ ${#completed_samples[@]} -gt 0 ]; then
        echo "partial_complete"
    else
        echo "none_complete"
    fi
}

# Function to create sample list
create_sample_list() {
    local fastq_dir=$1
    local sample_list_file=$2
    
    print_message $BLUE "Scanning for sample files in: $fastq_dir"
    
    # Find all _1.fastq.gz files and extract sample names
    local samples=($(find "$fastq_dir" -name "*_1.fastq.gz" -exec basename {} \; | sed 's/_1\.fastq\.gz$//' | sort))
    
    if [ ${#samples[@]} -eq 0 ]; then
        print_message $RED "Error: No sample files found in $fastq_dir"
        print_message $YELLOW "Expected format: SAMPLE_1.fastq.gz and SAMPLE_2.fastq.gz"
        return 1
    fi
    
    # Verify paired-end files exist for each sample
    local valid_samples=()
    for sample in "${samples[@]}"; do
        if [ -f "$fastq_dir/${sample}_1.fastq.gz" ] && [ -f "$fastq_dir/${sample}_2.fastq.gz" ]; then
            valid_samples+=("$sample")
            print_message $GREEN "✓ Valid sample: $sample"
        else
            print_message $RED "✗ Incomplete sample: $sample (missing paired-end file)"
        fi
    done
    
    if [ ${#valid_samples[@]} -eq 0 ]; then
        print_message $RED "Error: No valid paired-end samples found"
        return 1
    fi
    
    # Create sample list file
    printf "%s\n" "${valid_samples[@]}" > "$sample_list_file"
    
    print_message $GREEN "✓ Created sample list with ${#valid_samples[@]} samples: $sample_list_file"
    
    return ${#valid_samples[@]}
}

# Function to create filtered sample list (only incomplete samples)
create_filtered_sample_list() {
    local rdm_base_path=$1
    local original_sample_list=$2
    local filtered_sample_list=$3
    
    print_message $BLUE "Creating filtered sample list (incomplete samples only)..."
    
    local incomplete_samples=()
    
    # Check each sample and add only incomplete ones
    while IFS= read -r sample; do
        [ -z "$sample" ] && continue  # Skip empty lines
        
        local genotyped_vcf="${rdm_base_path}/5.Individual_VCF/${sample}_genotyped.vcf.gz"
        
        if [ ! -f "$genotyped_vcf" ] || [ ! -f "${genotyped_vcf}.tbi" ]; then
            incomplete_samples+=("$sample")
            print_message $YELLOW "  Adding to filtered list: $sample"
        else
            print_message $GREEN "  Skipping completed sample: $sample"
        fi
    done < "$original_sample_list"
    
    if [ ${#incomplete_samples[@]} -eq 0 ]; then
        print_message $GREEN "All samples are already complete!"
        return 1
    fi
    
    # Create filtered sample list
    printf "%s\n" "${incomplete_samples[@]}" > "$filtered_sample_list"
    
    print_message $GREEN "✓ Created filtered sample list with ${#incomplete_samples[@]} incomplete samples"
    
    return ${#incomplete_samples[@]}
}
