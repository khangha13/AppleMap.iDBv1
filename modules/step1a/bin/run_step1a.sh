#!/bin/bash
# =============================================================================
# STEP 1A MODULE - PER-SAMPLE VARIANT CALLING
# =============================================================================
# Called by bin/gatk_pipeline.sh

# Get script directory (module-specific to avoid clashes)
STEP1A_BIN_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STEP1A_MODULE_DIR="$(cd "${STEP1A_BIN_DIR}/.." && pwd)"
if [ -n "${PIPELINE_ROOT:-}" ]; then
    if [ -d "${PIPELINE_ROOT}" ]; then
        PIPELINE_ROOT="$(cd "${PIPELINE_ROOT}" && pwd)"
    else
        echo "[step1a] ⚠️  Provided PIPELINE_ROOT (${PIPELINE_ROOT}) does not exist; falling back to module-relative path." >&2
        PIPELINE_ROOT="$(cd "${STEP1A_BIN_DIR}/../../.." && pwd)"
    fi
else
    PIPELINE_ROOT="$(cd "${STEP1A_BIN_DIR}/../../.." && pwd)"
fi
export PIPELINE_ROOT

# Source required libraries
source "${PIPELINE_ROOT}/lib/logging.sh"
source "${PIPELINE_ROOT}/lib/slurm.sh"
source "${PIPELINE_ROOT}/lib/validation.sh"
source "${PIPELINE_ROOT}/lib/pipeline_common.sh"
source "${STEP1A_MODULE_DIR}/lib/functions.sh"
source "${PIPELINE_ROOT}/config/pipeline_config.sh"

# Guard: ensure pipeline_state_dir exists even if libraries weren’t loaded in sbatch spool
if ! command -v pipeline_state_dir >/dev/null 2>&1; then
    pipeline_state_dir() {
        local dataset_name="$1"
        if [ -n "${LOG_BASE_PATH:-}" ]; then
            echo "${LOG_BASE_PATH%/}/${dataset_name}"
        else
            echo "/tmp/${dataset_name}"
        fi
    }
fi

LAST_STEP1A_ARRAY_MAX=0

# Main Step 1A execution function
main() {
    local dataset_name="$1"
    local rdm_base_path="$2"
    shift 2

    # Optional overrides
    local custom_sample_name=""
    local custom_sample_list=""
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --sample)
                custom_sample_name="$2"; shift 2;;
            --sample-list)
                custom_sample_list="$2"; shift 2;;
            *)
                log_warn "Unknown option ignored: $1"
                shift;;
        esac
    done
    
    # Initialize logging
    init_logging "step1a" "pipeline" "${dataset_name}"
    
    log_info "Starting Step 1A execution"
    log_info "Dataset: $dataset_name"
    log_info "RDM Base Path: $rdm_base_path"
    
    # Validate inputs
    validate_step1a_inputs "$rdm_base_path"

    prepare_shared_reference_assets "${dataset_name}"
    
    local scratch_sample_dir="${PIPELINE_WORK_DIR%/}/step1a/${dataset_name}"
    mkdir -p "${scratch_sample_dir}"

    # Create or use sample list
    local sample_list_file=""
    if [ -n "${custom_sample_list}" ]; then
        if [ ! -f "${custom_sample_list}" ]; then
            log_error "Provided --sample-list not found: ${custom_sample_list}"
            exit 1
        fi
        sample_list_file="$(cd "$(dirname "${custom_sample_list}")" && pwd)/$(basename "${custom_sample_list}")"
        log_info "Using provided sample list: ${sample_list_file}"
    elif [ -n "${custom_sample_name}" ]; then
        # Validate FASTQ pair exists
        if [ ! -f "${rdm_base_path}/1.FASTQ/${custom_sample_name}_1.fastq.gz" ] || \
           [ ! -f "${rdm_base_path}/1.FASTQ/${custom_sample_name}_2.fastq.gz" ]; then
            log_error "FASTQ pair not found for sample '${custom_sample_name}' in ${rdm_base_path}/1.FASTQ"
            exit 1
        fi
        sample_list_file="${scratch_sample_dir}/sample_list_${dataset_name}_${custom_sample_name}_$(date +%Y%m%d_%H%M%S).txt"
        printf "%s\n" "${custom_sample_name}" > "${sample_list_file}"
        log_info "Created single-sample list for '${custom_sample_name}': ${sample_list_file}"
    else
        sample_list_file="${scratch_sample_dir}/sample_list_${dataset_name}_$(date +%Y%m%d_%H%M%S).txt"
        create_sample_list "$rdm_base_path/1.FASTQ" "$sample_list_file"
        # Do not rely on the function exit code (it returns the count); verify the file was created and non-empty
        if [ ! -s "$sample_list_file" ]; then
            log_error "Failed to create sample list (no samples found)"
            exit 1
        fi
    fi
    local sample_count
    sample_count=$(wc -l < "$sample_list_file")
    if [ "${sample_count}" -le 0 ]; then
        log_error "No valid samples detected in ${rdm_base_path}/1.FASTQ"
        exit 1
    fi
    
    # Check pipeline completion unless user specified an explicit list or sample
    if [ -z "${custom_sample_list}" ] && [ -z "${custom_sample_name}" ]; then
        local completion_status
        completion_status=$(check_pipeline_completion "$rdm_base_path" "$sample_list_file")
        case "$completion_status" in
            "all_complete")
                log_info "All samples are already complete"
                if confirm_action "Would you like to rerun Step 1A anyway?"; then
                    log_info "Proceeding with Step 1A rerun"
                else
                    log_info "Step 1A execution cancelled"
                    exit 0
                fi
                ;;
            "partial_complete")
                log_info "Some samples are already complete"
                if confirm_action "Would you like to run only incomplete samples?"; then
                    # Create filtered sample list
                    local filtered_sample_list="sample_list_${dataset_name}_incomplete_$(date +%Y%m%d_%H%M%S).txt"
                    create_filtered_sample_list "$rdm_base_path" "$sample_list_file" "$filtered_sample_list"
                    if [ $? -ne 0 ]; then
                        log_info "All samples are actually complete! Exiting."
                        exit 0
                    fi
                    sample_list_file="${filtered_sample_list}"
                    sample_count=$(wc -l < "$sample_list_file")
                    if [ "${sample_count}" -le 0 ]; then
                        log_info "No incomplete samples remain. Exiting."
                        exit 0
                    fi
                    log_info "Using filtered sample list with $sample_count incomplete samples"
                else
                    log_info "Proceeding with full Step 1A run"
                fi
                ;;
            "none_complete")
                log_info "No samples are complete - proceeding with full Step 1A run"
                ;;
        esac
    else
        log_info "Explicit sample selection provided; skipping completion prompts."
    fi
    
    # Create SLURM script
    # Normalise sample list path for the worker script
    sample_list_file="$(cd "$(dirname "${sample_list_file}")" && pwd)/$(basename "${sample_list_file}")"

    local slurm_script
    if ! slurm_script=$(create_step1a_slurm_script "$dataset_name" "$sample_count"); then
        log_error "Failed to generate Step 1A SLURM script"
        exit 1
    fi
    
    # Add delay to allow filesystem sync on network-mounted scratch directories
    log_info "Waiting 2 seconds for filesystem sync..."
    sleep 2
    

    # The file is guaranteed to exist if create_slurm_script returned successfully
    local array_range_max="${LAST_STEP1A_ARRAY_MAX:-$((sample_count - 1))}"
    
    # Submit job
    local job_id=""
    if job_id=$(submit_job "$slurm_script" "$rdm_base_path $sample_list_file" "$dataset_name" "1A"); then
        log_slurm_submission "$job_id" "$slurm_script" "$dataset_name" "$sample_count"

        local state_dir
        # Normalize state dir to LOG_BASE_PATH/<dataset> so master can find job IDs
        local state_dir="${LOG_BASE_PATH%/}/${dataset_name}"
        mkdir -p "${state_dir}"
        printf '%s\n' "${job_id}" > "${state_dir}/step1a_job_id.txt"
        printf '%s\n' "${sample_list_file}" > "${state_dir}/step1a_samples.list"
        rm -f "${state_dir}/step1a_failed.flag" 2>/dev/null || true
        touch "${state_dir}/step1a_running.flag"
        log_info "Recorded Step 1A job ID to ${state_dir}/step1a_job_id.txt"

        echo "✓ Step 1A submitted successfully!"
        echo "Job ID: $job_id"
        echo "Array jobs: 0-${array_range_max}"
        echo "Monitor with: squeue -u \$USER"
        echo "Cancel with: scancel $job_id"
        echo "View logs: ls ${LOG_BASE_PATH}/${dataset_name}/"
    else
        log_error "Failed to submit Step 1A"
        exit 1
    fi
}

# Function to execute Step 1A pipeline for a single sample
execute_step1a_pipeline() {
    local sample="$1"
    local rdm_base_path="$2"
    local dataset_name="${3:-$(basename "${rdm_base_path}")}"
    
    log_info "Executing Step 1A pipeline for sample: $sample"
    
    # Set up environment variables
    local NUM_THREADS="${STEP1A_CPUS_PER_TASK:-10}"
    local MEMORY="${STEP1A_MEMORY:-32G}"
    # Get reference paths (don't declare as local to avoid shadowing global REFERENCE_GENOME)
    local ref_genome known_sites adapter_file
    ref_genome="$(get_reference_fasta)"
    known_sites="$(get_known_sites_vcf)"
    adapter_file="$(get_adapter_fasta)"
    
    # Use local variables for the rest of the function
    local REFERENCE_GENOME="$ref_genome"
    local KNOWN_SITES="$known_sites"
    local ADAPTER_FILE="$adapter_file"
    local READ_GROUP
    READ_GROUP="$(build_read_group "$sample")"

    if [ ! -f "$REFERENCE_GENOME" ]; then
        error_exit "Reference genome not found at $REFERENCE_GENOME"
    fi
    if [ ! -f "$KNOWN_SITES" ]; then
        error_exit "Known sites VCF not found at $KNOWN_SITES"
    fi
    if [ ! -f "$ADAPTER_FILE" ]; then
        error_exit "Adapter FASTA not found at $ADAPTER_FILE"
    fi
    if [ ! -f "${TRIMMOMATIC_JAR}" ]; then
        error_exit "Trimmomatic jar not found at ${TRIMMOMATIC_JAR}"
    fi
    
    ensure_sample_directories "$rdm_base_path" "$sample"
    
    # Determine resume step
    local completed_step
    completed_step=$(detect_completed_step "$sample" "$rdm_base_path")
    if [ "$completed_step" -ge 7 ]; then
        log_info "All steps already completed for ${sample}; skipping execution."
        return 0
    fi
    if [ "$completed_step" -gt 0 ]; then
        log_warn "Existing outputs detected for ${sample} (last completed step: ${completed_step}). Steps will be rerun to ensure consistency."
    fi
    STARTING_STEP=1
    log_info "Starting pipeline from step ${STARTING_STEP}."
    
    # Change to working directory
    local CLEANUP_TMPDIR="false"
    local WORK_TMPDIR
    if [ -n "${TMPDIR:-}" ]; then
        WORK_TMPDIR="${TMPDIR}"
    else
        WORK_TMPDIR="$(mktemp -d "${SCRATCH_BASE_PATH%/}/step1a_${sample}_XXXXXX")"
        CLEANUP_TMPDIR="true"
    fi
    mkdir -p "${WORK_TMPDIR}"
    cd "${WORK_TMPDIR}" || error_exit "Failed to change to working directory (${WORK_TMPDIR})"
    
    # Restore from backup if starting from step >= 2 and backup exists (allows Step 2 restore)
    if [ "$STARTING_STEP" -ge 2 ]; then
        if restore_from_backup "$sample" "$dataset_name"; then
            log_info "Restored files from backup, continuing from step $STARTING_STEP"
            # Infer local starting step based on presence of trimmed outputs
            local inferred_step
            inferred_step="$(infer_starting_step_from_local "$sample" 2>/dev/null || echo 1)"
            if [ "$inferred_step" -ge 3 ]; then
                STARTING_STEP="$inferred_step"
                log_info "Detected restored trimmed outputs; starting from Step ${STARTING_STEP}."
            fi
        else
            log_info "No backup found, will proceed with normal file copying"
        fi
    fi
    
    # Copy all reference files to working directory (matching original script approach)
    # This ensures all files needed for calculation are local and accessible
    # Uses symlinks if shared reference files exist (for array jobs)
    local ref_files_result
    ref_files_result="$(copy_reference_files_to_workdir "${REFERENCE_GENOME}" "${KNOWN_SITES}" "${ADAPTER_FILE}" "${dataset_name}")"
    
    # Parse returned basenames (format: ref_basename|known_sites_basename|adapter_basename)
    local LOCAL_REF_GENOME LOCAL_KNOWN_SITES LOCAL_ADAPTER_FILE
    IFS='|' read -r LOCAL_REF_GENOME LOCAL_KNOWN_SITES LOCAL_ADAPTER_FILE <<< "${ref_files_result}"
    
    # Copy input FASTQ files using rsync (matching original script approach)
    # Only if not restored from backup or if FASTQ files are missing
    if [ ! -f "${sample}_forward_paired.fastq.gz" ] && [ ! -f "${sample}_reverse_paired.fastq.gz" ]; then
        if [ ! -f "${rdm_base_path}/1.FASTQ/${sample}_1.fastq.gz" ] || [ ! -f "${rdm_base_path}/1.FASTQ/${sample}_2.fastq.gz" ]; then
            error_exit "Missing FASTQ pair for sample ${sample} in ${rdm_base_path}/1.FASTQ"
        fi
        rsync -rhivPt "${rdm_base_path}/1.FASTQ/${sample}_1.fastq.gz" "${rdm_base_path}/1.FASTQ/${sample}_2.fastq.gz" . || error_exit "Failed to copy FASTQ files to working directory"
    else
        log_info "FASTQ files already present (restored from backup or previous step)"
    fi
    
    # Step 1: FastQC on raw reads
    if verify_starting_step 1; then
        start_step_timer "Step 1: FastQC on raw reads"
        run_fastqc "$sample" "$NUM_THREADS"
        end_step_timer "Step 1: FastQC on raw reads"
        
        # Copy results to RDM
        rsync -rhivPt *.html "${rdm_base_path}/2.FASTQC_pre_trimmed/${sample}/" || error_exit "Failed to copy FastQC HTML results to RDM"
        rsync -rhivPt *.zip "${rdm_base_path}/2.FASTQC_pre_trimmed/${sample}/" || error_exit "Failed to copy FastQC ZIP results to RDM"
        rm *.html *.zip
    fi
    
    # Step 2: Read trimming with Trimmomatic
    if verify_starting_step 2; then
        start_step_timer "Step 2: Read trimming with Trimmomatic"
        run_trimmomatic "$sample" "$MEMORY" "$NUM_THREADS" "$LOCAL_ADAPTER_FILE"
        end_step_timer "Step 2: Read trimming with Trimmomatic"
        
        # Step 2b: FastQC on trimmed reads
        start_step_timer "Step 2b: FastQC on trimmed reads"
        run_fastqc "$sample" "$NUM_THREADS" "${sample}_forward_paired.fastq.gz" "${sample}_reverse_paired.fastq.gz"
        end_step_timer "Step 2b: FastQC on trimmed reads"
        
        # Copy results to RDM
        rsync -rhivPt *.html "${rdm_base_path}/3.FASTQC_post_trimmed/${sample}/" || error_exit "Failed to copy FastQC HTML results to RDM"
        rsync -rhivPt *.zip "${rdm_base_path}/3.FASTQC_post_trimmed/${sample}/" || error_exit "Failed to copy FastQC ZIP results to RDM"
        rm *.html *.zip
    fi
    
    # Step 3: BWA alignment
    # If resuming at Step >=3, ensure post-trim FastQC exists on RDM; regenerate if missing
    if [ "$STARTING_STEP" -ge 3 ]; then
        local post_trim_dir="${rdm_base_path}/3.FASTQC_post_trimmed/${sample}"
        if [ ! -d "${post_trim_dir}" ] || ! ls "${post_trim_dir}"/*.html >/dev/null 2>&1; then
            start_step_timer "Step 2b: FastQC on trimmed reads (resume)"
            run_fastqc "$sample" "$NUM_THREADS" "${sample}_forward_paired.fastq.gz" "${sample}_reverse_paired.fastq.gz"
            end_step_timer "Step 2b: FastQC on trimmed reads (resume)"
            rsync -rhivPt *.html "${post_trim_dir}/" || error_exit "Failed to copy FastQC HTML results (resume) to RDM"
            rsync -rhivPt *.zip "${post_trim_dir}/" || error_exit "Failed to copy FastQC ZIP results (resume) to RDM"
            rm -f *.html *.zip
        fi
    fi
    # Backup trimmed outputs after Step 2b (paired FASTQs only)
    backup_current_state "Step 2: Trimmomatic" "$sample" "$dataset_name"

    if verify_starting_step 3; then
        start_step_timer "Step 3: BWA alignment"
        run_bwa_alignment "$sample" "$NUM_THREADS" "$LOCAL_REF_GENOME" "$READ_GROUP"
        end_step_timer "Step 3: BWA alignment"
        
        # Step 3b: BAM indexing
        start_step_timer "Step 3b: BAM indexing"
        "${SAMTOOLS_PATH}" index "${sample}_sorted.bam"
        end_step_timer "Step 3b: BAM indexing"
        
        # Backup after Step 3
        backup_current_state "Step 3: BWA alignment" "$sample" "$dataset_name"
    fi
    
    # Step 4: Duplicate marking
    if verify_starting_step 4; then
        start_step_timer "Step 4: Duplicate marking"
        run_mark_duplicates "$sample" "$MEMORY"
        end_step_timer "Step 4: Duplicate marking"
        
        # Step 4b: Deduplicated BAM indexing
        start_step_timer "Step 4b: Deduplicated BAM indexing"
        "${SAMTOOLS_PATH}" index "${sample}_dedup.bam"
        end_step_timer "Step 4b: Deduplicated BAM indexing"
        
        # Backup after Step 4
        backup_current_state "Step 4: Duplicate marking" "$sample" "$dataset_name"
    fi
    
    # Step 5: Base quality score recalibration
    if verify_starting_step 5; then
        start_step_timer "Step 5a: Base recalibration table generation"
        run_base_recalibration "$sample" "$MEMORY" "$LOCAL_REF_GENOME" "$LOCAL_KNOWN_SITES"
        end_step_timer "Step 5a: Base recalibration table generation"
        
        start_step_timer "Step 5b: Apply base recalibration"
        run_apply_bqsr "$sample" "$MEMORY" "$LOCAL_REF_GENOME"
        end_step_timer "Step 5b: Apply base recalibration"
        
        # Index recalibrated BAM
        "${SAMTOOLS_PATH}" index "${sample}_recal.bam"
        
        # Backup after Step 5
        backup_current_state "Step 5: Base quality score recalibration" "$sample" "$dataset_name"
    fi
    
    # Step 6: Variant calling with HaplotypeCaller
    if verify_starting_step 6; then
        start_step_timer "Step 6: Variant calling with HaplotypeCaller"
        run_haplotype_caller "$sample" "$MEMORY" "$LOCAL_REF_GENOME" "$NUM_THREADS"
        end_step_timer "Step 6: Variant calling with HaplotypeCaller"
        
        # Backup after Step 6
        backup_current_state "Step 6: Variant calling with HaplotypeCaller" "$sample" "$dataset_name"
    fi
    
    # Step 7: Genotype GVCF to VCF conversion
    if verify_starting_step 7; then
        start_step_timer "Step 7: Genotype GVCF to VCF conversion"
        run_genotype_gvcfs "$sample" "$MEMORY" "$LOCAL_REF_GENOME"
        end_step_timer "Step 7: Genotype GVCF to VCF conversion"
    fi
    
    # Copy results to RDM
    copy_results_to_rdm "$sample" "$rdm_base_path"
    
    log_info "Step 1A pipeline completed successfully for sample: $sample"

    if [ "${CLEANUP_TMPDIR}" = "true" ]; then
        rm -rf "${WORK_TMPDIR}"
    fi
}

# Validate Step 1A inputs
validate_step1a_inputs() {
    local rdm_base_path="$1"
    
    log_info "Validating Step 1A inputs"
    
    # Check if RDM base path exists
    if [ ! -d "$rdm_base_path" ]; then
        log_error "RDM base path does not exist: $rdm_base_path"
        exit 1
    fi
    
    # Check if 1.FASTQ directory exists
    if [ ! -d "$rdm_base_path/1.FASTQ" ]; then
        log_error "1.FASTQ directory does not exist: $rdm_base_path/1.FASTQ"
        exit 1
    fi
    
    # Validate directory structure
    validate_directory_structure "$rdm_base_path"
    
    log_info "Step 1A input validation completed successfully"
}

# Create Step 1A SLURM script
create_step1a_slurm_script() {
    local dataset_name="$1"
    local sample_count="$2"
    local array_max=$((sample_count - 1))
    
    log_info "Creating Step 1A SLURM script"
    
    # Get Step 1A configuration
    local config
    config=$(get_step1a_config)
    declare -A config_map=()
    while IFS='=' read -r key value; do
        [ -z "$key" ] && continue
        config_map["$key"]="$value"
    done <<< "$config"

    local configured_array_limit="${config_map[ARRAY_MAX]:-}"
    if [[ -n "$configured_array_limit" && "$configured_array_limit" -gt 0 && "$array_max" -ge "$configured_array_limit" ]]; then
        log_warn "Sample count (${sample_count}) exceeds configured array limit (${configured_array_limit}); truncating array range."
        array_max=$((configured_array_limit - 1))
    fi
    LAST_STEP1A_ARRAY_MAX="$array_max"
    
    # Create SLURM script
    mkdir -p "${PIPELINE_SLURM_SCRIPT_DIR}"
    local slurm_script="${PIPELINE_SLURM_SCRIPT_DIR}/Apple_GATK_1A_${dataset_name}_$(date +%Y%m%d_%H%M%S).sh"
    
    # Create configuration string
    local config_string="job_name=Apple_GATK_1A_${dataset_name}
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
    
    # Create SLURM script using template
    local template="${STEP1A_MODULE_DIR}/templates/step1a_array.sh"
    if ! create_slurm_script "${template}" "$config_string" "$slurm_script" "$dataset_name" "1A" "${PIPELINE_ROOT}"; then
        log_error "Failed to create SLURM script"
        return 1
    fi
    
    echo "$slurm_script"
}

# Run main function when executed directly
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    main "$@"
fi
