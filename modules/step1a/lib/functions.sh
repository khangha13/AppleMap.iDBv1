#!/bin/bash
# =============================================================================
# STEP 1A FUNCTIONS
# =============================================================================
# Step 1A specific functions for per-sample variant calling
# Extracted from Apple_GATK_Pipeline_1a_CallVariantsPerSampleKHv2.sh

# Utility: ensure per-sample directories exist on RDM
ensure_sample_directories() {
    local rdm_base="$1"
    local sample="$2"

    mkdir -p \
        "${rdm_base}/2.FASTQC_pre_trimmed/${sample}" \
        "${rdm_base}/3.FASTQC_post_trimmed/${sample}" \
        "${rdm_base}/4.BAM/${sample}" \
        "${rdm_base}/5.Individual_VCF" \
        "${rdm_base}/9.Metrics/${sample}"
}

# Utility: setup shared reference files for array jobs (task 0 only)
setup_shared_reference_files() {
    local dataset_name="$1"
    local reference_genome="$2"
    local known_sites="$3"
    local adapter_file="$4"
    
    # Only run if this is task 0
    if [ "${SLURM_ARRAY_TASK_ID:-}" != "0" ]; then
        return 0
    fi
    
    log_info "Setting up shared reference files for all array jobs"
    
    # Create shared directories
    local shared_base="${SCRATCH_BASE_PATH%/}/${dataset_name}_shared"
    local shared_ref_dir="${shared_base}/Reference_genome"
    local shared_known_dir="${shared_base}/Known_sites"
    local shared_adapter_dir="${shared_base}/Adapter_file"
    
    mkdir -p "${shared_ref_dir}" "${shared_known_dir}" "${shared_adapter_dir}"
    
    # Get basenames
    local ref_basename known_sites_basename adapter_basename
    ref_basename="$(basename "${reference_genome}")"
    known_sites_basename="$(basename "${known_sites}")"
    adapter_basename="$(basename "${adapter_file}")"
    
    # Check if reference files already exist in shared location
    if [ -f "${shared_ref_dir}/${ref_basename}" ]; then
        log_info "Shared reference files already exist, skipping copy"
    else
        log_info "Copying reference genome and indices to shared directory"
        rsync -rhivPt "${reference_genome}" "${shared_ref_dir}/" || error_exit "Failed to copy reference genome to shared directory"
        
        # Copy reference index (.fai) if it exists
        if [ -f "${reference_genome}.fai" ]; then
            rsync -rhivPt "${reference_genome}.fai" "${shared_ref_dir}/" || log_info "Warning: Failed to copy reference genome index (.fai), continuing anyway"
        fi
        
        # Copy BWA index files if they exist
        for ext in amb ann bwt pac sa alt; do
            if [ -f "${reference_genome}.${ext}" ]; then
                rsync -rhivPt "${reference_genome}.${ext}" "${shared_ref_dir}/" || log_info "Warning: Failed to copy reference genome BWA index (.${ext}), continuing anyway"
            fi
        done
        
        # Copy reference dictionary (.dict) if it exists
        local ref_dict="${reference_genome%.*}.dict"
        if [ -f "${ref_dict}" ]; then
            rsync -rhivPt "${ref_dict}" "${shared_ref_dir}/" || log_info "Warning: Failed to copy reference genome dictionary (.dict), continuing anyway"
        fi
        
        log_info "Copying known sites to shared directory"
        rsync -rhivPt "${known_sites}" "${shared_known_dir}/" || error_exit "Failed to copy known sites VCF to shared directory"
        
        # Copy known sites index (.tbi) if it exists
        if [ -f "${known_sites}.tbi" ]; then
            rsync -rhivPt "${known_sites}.tbi" "${shared_known_dir}/" || log_info "Warning: Failed to copy known sites VCF index (.tbi), continuing anyway"
        fi
        # Copy known sites index (.idx) if it exists
        if [ -f "${known_sites}.idx" ]; then
            rsync -rhivPt "${known_sites}.idx" "${shared_known_dir}/" || log_info "Warning: Failed to copy known sites VCF index (.idx), continuing anyway"
        fi
        
        log_info "Copying adapter file to shared directory"
        rsync -rhivPt "${adapter_file}" "${shared_adapter_dir}/" || error_exit "Failed to copy adapter file to shared directory"
        
        # Verify all files exist before creating marker
        if [ -f "${shared_ref_dir}/${ref_basename}" ] && \
           [ -f "${shared_known_dir}/${known_sites_basename}" ] && \
           [ -f "${shared_adapter_dir}/${adapter_basename}" ]; then
            # Create marker file AFTER all files are verified to exist
            touch "${shared_base}/.initialized"
            log_info "Shared reference files setup completed"
        else
            error_exit "Shared reference files setup incomplete - some files missing"
        fi
    fi
    
    # Ensure known sites indexes exist in shared directory even if files already existed previously
    if [ -f "${shared_known_dir}/${known_sites_basename}" ]; then
        if [ ! -f "${shared_known_dir}/${known_sites_basename}.tbi" ] && [ -f "${known_sites}.tbi" ]; then
            rsync -rhivPt "${known_sites}.tbi" "${shared_known_dir}/" || log_info "Warning: Failed to copy known sites VCF index (.tbi), continuing anyway"
        fi
        if [ ! -f "${shared_known_dir}/${known_sites_basename}.idx" ] && [ -f "${known_sites}.idx" ]; then
            rsync -rhivPt "${known_sites}.idx" "${shared_known_dir}/" || log_info "Warning: Failed to copy known sites VCF index (.idx), continuing anyway"
        fi
    fi
}

# Utility: get paths to shared reference files if they exist
get_shared_reference_paths() {
    local dataset_name="$1"
    local reference_genome="$2"
    local known_sites="$3"
    local adapter_file="$4"
    
    # Only check for shared files in array job context
    if [ -z "${SLURM_ARRAY_TASK_ID:-}" ]; then
        echo ""
        return 0
    fi
    
    local shared_base="${SCRATCH_BASE_PATH%/}/${dataset_name}_shared"
    local shared_ref_dir="${shared_base}/Reference_genome"
    local shared_known_dir="${shared_base}/Known_sites"
    local shared_adapter_dir="${shared_base}/Adapter_file"
    
    local ref_basename known_sites_basename adapter_basename
    ref_basename="$(basename "${reference_genome}")"
    known_sites_basename="$(basename "${known_sites}")"
    adapter_basename="$(basename "${adapter_file}")"
    
    # Check if shared files exist
    if [ -f "${shared_ref_dir}/${ref_basename}" ] && \
       [ -f "${shared_known_dir}/${known_sites_basename}" ] && \
       [ -f "${shared_adapter_dir}/${adapter_basename}" ]; then
        echo "${shared_ref_dir}/${ref_basename}|${shared_known_dir}/${known_sites_basename}|${shared_adapter_dir}/${adapter_basename}"
    else
        echo ""
    fi
}

# Utility: wait for shared reference files to be initialized (for non-task-0 jobs)
wait_for_shared_reference_files() {
    local dataset_name="$1"
    local reference_genome="$2"
    local known_sites="$3"
    local adapter_file="$4"
    
    # Only wait if this is not task 0
    if [ "${SLURM_ARRAY_TASK_ID:-}" = "0" ]; then
        return 0
    fi
    
    local shared_base="${SCRATCH_BASE_PATH%/}/${dataset_name}_shared"
    local marker_file="${shared_base}/.initialized"
    local timeout="${PIPELINE_SHARED_REF_TIMEOUT:-120}"
    local poll_interval="${PIPELINE_SHARED_REF_POLL_INTERVAL:-3}"
    local elapsed=0
    
    log_info "Waiting for shared reference files initialization (timeout: ${timeout}s, polling every ${poll_interval}s)"
    
    while [ $elapsed -lt $timeout ]; do
        # Check marker file exists
        if [ -f "${marker_file}" ]; then
            # Verify actual reference files exist (not just marker)
            local ref_basename known_sites_basename adapter_basename
            ref_basename="$(basename "${reference_genome}")"
            known_sites_basename="$(basename "${known_sites}")"
            adapter_basename="$(basename "${adapter_file}")"
            
            if [ -f "${shared_base}/Reference_genome/${ref_basename}" ] && \
               [ -f "${shared_base}/Known_sites/${known_sites_basename}" ] && \
               [ -f "${shared_base}/Adapter_file/${adapter_basename}" ]; then
                log_info "Shared reference files are ready (waited ${elapsed}s)"
                return 0
            fi
        fi
        
        sleep "$poll_interval"
        elapsed=$((elapsed + poll_interval))
        
        # Log progress every 15 seconds
        if [ $((elapsed % 15)) -eq 0 ]; then
            log_info "Still waiting for shared reference files... (${elapsed}/${timeout}s)"
        fi
    done
    
    error_exit "Timeout waiting for shared reference files initialization (${timeout}s exceeded)"
}

# Utility: copy all reference files to working directory (matching original script approach)
# Uses symlinks if shared reference files exist, otherwise copies directly
copy_reference_files_to_workdir() {
    local reference_genome="$1"
    local known_sites="$2"
    local adapter_file="$3"
    local dataset_name="${4:-}"
    
    # Check if shared reference files exist (for array jobs)
    local shared_paths
    if [ -n "${dataset_name}" ] && [ -n "${SLURM_ARRAY_TASK_ID:-}" ]; then
        shared_paths="$(get_shared_reference_paths "${dataset_name}" "${reference_genome}" "${known_sites}" "${adapter_file}")"
    else
        shared_paths=""
    fi
    
    local ref_basename known_sites_basename adapter_basename
    
    if [ -n "${shared_paths}" ]; then
        # Use symlinks for shared reference files
        log_info "Using shared reference files (creating symlinks)"
        
        IFS='|' read -r shared_ref shared_known shared_adapter <<< "${shared_paths}"
        
        ref_basename="$(basename "${shared_ref}")"
        known_sites_basename="$(basename "${shared_known}")"
        adapter_basename="$(basename "${shared_adapter}")"
        
        # Create symlinks
        ln -sf "${shared_ref}" "${ref_basename}" || error_exit "Failed to create symlink for reference genome"
        
        # Create symlinks for index files if they exist
        if [ -f "${shared_ref}.fai" ]; then
            ln -sf "${shared_ref}.fai" "${ref_basename}.fai" || log_info "Warning: Failed to create symlink for reference genome index (.fai), continuing anyway"
        fi
        # Create symlinks for BWA index files if they exist
        for ext in amb ann bwt pac sa alt; do
            if [ -f "${shared_ref}.${ext}" ]; then
                ln -sf "${shared_ref}.${ext}" "${ref_basename}.${ext}" || log_info "Warning: Failed to create symlink for BWA index (.${ext}), continuing anyway"
            fi
        done
        
        local ref_dict="${shared_ref%.*}.dict"
        if [ -f "${ref_dict}" ]; then
            local ref_dict_basename
            ref_dict_basename="$(basename "${ref_dict}")"
            ln -sf "${ref_dict}" "${ref_dict_basename}" || log_info "Warning: Failed to create symlink for reference genome dictionary (.dict), continuing anyway"
        fi
        
        ln -sf "${shared_known}" "${known_sites_basename}" || error_exit "Failed to create symlink for known sites VCF"
        
        if [ -f "${shared_known}.tbi" ]; then
            ln -sf "${shared_known}.tbi" "${known_sites_basename}.tbi" || log_info "Warning: Failed to create symlink for known sites VCF index (.tbi), continuing anyway"
        fi
        if [ -f "${shared_known}.idx" ]; then
            ln -sf "${shared_known}.idx" "${known_sites_basename}.idx" || log_info "Warning: Failed to create symlink for known sites VCF index (.idx), continuing anyway"
        fi
        
        ln -sf "${shared_adapter}" "${adapter_basename}" || error_exit "Failed to create symlink for adapter file"
        
        log_info "Reference files linked successfully"
    else
        # Copy files directly (non-array job or shared files not available)
        log_info "Copying reference files to working directory"
        
        ref_basename="$(basename "${reference_genome}")"
        rsync -rhivPt "${reference_genome}" "${ref_basename}" || error_exit "Failed to copy reference genome to working directory"
        
        # Copy reference index (.fai) if it exists
        if [ -f "${reference_genome}.fai" ]; then
            rsync -rhivPt "${reference_genome}.fai" "${ref_basename}.fai" || log_info "Warning: Failed to copy reference genome index (.fai), continuing anyway"
        fi
        
        # Copy reference dictionary (.dict) if it exists
        local ref_dict="${reference_genome%.*}.dict"
        if [ -f "${ref_dict}" ]; then
            local ref_dict_basename
            ref_dict_basename="$(basename "${ref_dict}")"
            rsync -rhivPt "${ref_dict}" "${ref_dict_basename}" || log_info "Warning: Failed to copy reference genome dictionary (.dict), continuing anyway"
        fi
        
        # Copy BWA index files if present
        for ext in amb ann bwt pac sa alt; do
            if [ -f "${reference_genome}.${ext}" ]; then
                rsync -rhivPt "${reference_genome}.${ext}" "${ref_basename}.${ext}" || log_info "Warning: Failed to copy BWA index (.${ext}), continuing anyway"
            fi
        done
        
        # Copy known sites VCF and its index
        known_sites_basename="$(basename "${known_sites}")"
        rsync -rhivPt "${known_sites}" "${known_sites_basename}" || error_exit "Failed to copy known sites VCF to working directory"
        
        # Copy known sites index (.tbi) if it exists
        if [ -f "${known_sites}.tbi" ]; then
            rsync -rhivPt "${known_sites}.tbi" "${known_sites_basename}.tbi" || log_info "Warning: Failed to copy known sites VCF index (.tbi), continuing anyway"
        fi
        # Copy known sites index (.idx) if it exists
        if [ -f "${known_sites}.idx" ]; then
            rsync -rhivPt "${known_sites}.idx" "${known_sites_basename}.idx" || log_info "Warning: Failed to copy known sites VCF index (.idx), continuing anyway"
        fi
        
        # Copy adapter file
        adapter_basename="$(basename "${adapter_file}")"
        rsync -rhivPt "${adapter_file}" "${adapter_basename}" || error_exit "Failed to copy adapter file to working directory"
        
        log_info "Reference files copied successfully"
    fi
    
    # Return basenames for use in pipeline functions
    echo "${ref_basename}|${known_sites_basename}|${adapter_basename}"
}

# Infer local starting step from presence of trimmed outputs in the current working directory
infer_starting_step_from_local() {
    local sample="$1"
    if [ -s "${sample}_forward_paired.fastq.gz" ] && [ -s "${sample}_reverse_paired.fastq.gz" ]; then
        echo 3
        return 0
    fi
    echo 1
}

# Utility: construct read group string
build_read_group() {
    local sample="$1"
    printf '%s\n' "@RG\\tID:${sample}\\tSM:${sample}\\tPL:ILLUMINA\\tLB:LIB1\\tPU:${sample}"
}

# Utility: determine the highest completed step (1-7) based on RDM artefacts
detect_completed_step() {
    local sample="$1"
    local rdm_base="$2"

    local vcf_dir="${rdm_base}/5.Individual_VCF"
    local bam_dir="${rdm_base}/4.BAM/${sample}"
    local qc_post_dir="${rdm_base}/3.FASTQC_post_trimmed/${sample}"
    local qc_pre_dir="${rdm_base}/2.FASTQC_pre_trimmed/${sample}"

    if [ -f "${vcf_dir}/${sample}_genotyped.vcf.gz" ]; then
        echo 7
        return
    fi
    if [ -f "${vcf_dir}/${sample}_raw.g.vcf.gz" ]; then
        echo 6
        return
    fi
    if [ -f "${bam_dir}/${sample}_recal.bam" ]; then
        echo 5
        return
    fi
    if [ -f "${bam_dir}/${sample}_dedup.bam" ]; then
        echo 4
        return
    fi
    if [ -f "${bam_dir}/${sample}_sorted.bam" ]; then
        echo 3
        return
    fi
    if ls "${qc_post_dir}/"* 1>/dev/null 2>&1; then
        echo 2
        return
    fi
    if ls "${qc_pre_dir}/"* 1>/dev/null 2>&1; then
        echo 1
        return
    fi
    echo 0
}

# Function to start timing a step with enhanced formatting
start_step_timer() {
    local step_name="$1"
    STEP_START_TIME=$(date +%s)
    CURRENT_STEP_NAME="$step_name"
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    
    # Log step start with professional logging
    log_info "Starting step: $step_name"
    
    echo ""
    echo "╔══════════════════════════════════════════════════════════════════════════════════╗"
    echo "║                           STARTING: $step_name"
    echo "║                           Started at: $timestamp"
    echo "╚══════════════════════════════════════════════════════════════════════════════════╝"
    echo ""
    
    # Log step start to file
    log_info "Step started: $step_name at $timestamp"
}

# Function to end timing a step and log duration with enhanced formatting
end_step_timer() {
    local step_name="$1"
    STEP_END_TIME=$(date +%s)
    local duration=$((STEP_END_TIME - STEP_START_TIME))
    local hours=$((duration / 3600))
    local minutes=$(((duration % 3600) / 60))
    local seconds=$((duration % 60))
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    
    # Log step completion with professional logging
    if [ $hours -gt 0 ]; then
        log_info "Completed step: $step_name (Duration: ${hours}h ${minutes}m ${seconds}s)"
    elif [ $minutes -gt 0 ]; then
        log_info "Completed step: $step_name (Duration: ${minutes}m ${seconds}s)"
    else
        log_info "Completed step: $step_name (Duration: ${seconds}s)"
    fi
    
    echo ""
    echo "╔══════════════════════════════════════════════════════════════════════════════════╗"
    echo "║                           COMPLETED: $step_name"
    echo "║                           Completed at: $timestamp"
    
    if [ $hours -gt 0 ]; then
        echo "║                           Duration: ${hours}h ${minutes}m ${seconds}s"
    elif [ $minutes -gt 0 ]; then
        echo "║                           Duration: ${minutes}m ${seconds}s"
    else
        echo "║                           Duration: ${seconds}s"
    fi
    echo "╚══════════════════════════════════════════════════════════════════════════════════╝"
    echo ""
}

# Function to verify starting step
verify_starting_step() {
    local step_number="$1"
    log_message "Verifying starting step: $step_number"
    
    if [ "$STARTING_STEP" -gt "$step_number" ]; then
        log_message "Step $step_number: SKIPPED (already completed)"
        return 1
    else
        log_message "Step $step_number: WILL RUN"
        return 0
    fi
}

# Function to run FastQC on provided FASTQ files (or all *.fastq.gz if none specified)
run_fastqc() {
    local sample="$1"
    local num_threads="$2"
    shift 2
    local fastq_inputs=("$@")
    
    if [ ${#fastq_inputs[@]} -eq 0 ]; then
        fastq_inputs=(*.fastq.gz)
    fi
    if [ ${#fastq_inputs[@]} -eq 0 ]; then
        error_exit "No FASTQ files found for FastQC in working directory"
    fi
    
    log_message "Running FastQC for sample: $sample on ${#fastq_inputs[@]} files"
    
    if ! "${FASTQC_PATH}" -t "$num_threads" "${fastq_inputs[@]}"; then
        error_exit "FastQC failed"
    fi
    
    log_message "FastQC completed successfully"
}

# Function to run Trimmomatic
run_trimmomatic() {
    local sample="$1"
    local memory="$2"
    local num_threads="$3"
    local adapter_file="$4"
    
    log_message "Running Trimmomatic for sample: $sample"
    
    # Adapter file should already be in working directory (copied by copy_reference_files_to_workdir)
    # Use basename for ILLUMINACLIP parameter
    local adapter_basename
    adapter_basename="$(basename "${adapter_file}")"
    
    # Build ILLUMINACLIP parameter: adapter_file:settings
    # Format: ILLUMINACLIP:adapter_file:mismatches:palindrome:simple:minAdapterLength:keepBothReads
    # Matching original: ILLUMINACLIP:$(basename "$ADAPTER_FILE"):${ILLUMINA_CLIP_SETTINGS}
    local illuminaclip_param="ILLUMINACLIP:${adapter_basename}:${TRIMMOMATIC_ILLUMINACLIP_SETTINGS}"
    
    if ! java -Xmx"${memory}" -jar "${TRIMMOMATIC_JAR}" PE \
        -threads "$num_threads" \
        *1.fastq.gz *2.fastq.gz \
        "${sample}_forward_paired.fastq.gz" "${sample}_forward_unpaired.fastq.gz" \
        "${sample}_reverse_paired.fastq.gz" "${sample}_reverse_unpaired.fastq.gz" \
        "${illuminaclip_param}" \
        LEADING:"${TRIMMOMATIC_LEADING}" \
        TRAILING:"${TRIMMOMATIC_TRAILING}" \
        SLIDINGWINDOW:"${TRIMMOMATIC_SLIDINGWINDOW}" \
        AVGQUAL:"${TRIMMOMATIC_AVGQUAL}" \
        MINLEN:"${TRIMMOMATIC_MINLEN}"; then
        error_exit "Trimmomatic failed"
    fi
    
    log_message "Trimmomatic completed successfully"
}

# Function to run BWA alignment
run_bwa_alignment() {
    local sample="$1"
    local num_threads="$2"
    local reference_genome="$3"
    local read_group="$4"
    
    log_message "Running BWA alignment for sample: $sample"
    
    if ! "${BWA_PATH}" mem ${BWA_MEM_FLAGS} -t "$num_threads" -R "$read_group" \
        "$reference_genome" \
        "${sample}_forward_paired.fastq.gz" \
        "${sample}_reverse_paired.fastq.gz" \
        | "${SAMTOOLS_PATH}" sort -o "${sample}_sorted.bam"; then
        error_exit "BWA alignment failed"
    fi
    
    log_message "BWA alignment completed successfully"
}

# Function to run GATK MarkDuplicates
run_mark_duplicates() {
    local sample="$1"
    local memory="$2"
    
    log_message "Running GATK MarkDuplicates for sample: $sample"
    
    if ! "${GATK_COMMAND}" --java-options "-Xmx${memory}" MarkDuplicates \
        -I "${sample}_sorted.bam" \
        -O "${sample}_dedup.bam" \
        -M "${sample}_dedup_metrics_file.txt"; then
        error_exit "GATK MarkDuplicates failed"
    fi
    
    log_message "GATK MarkDuplicates completed successfully"
}

# Function to run GATK BaseRecalibrator
run_base_recalibration() {
    local sample="$1"
    local memory="$2"
    local reference_genome="$3"
    local known_sites="$4"
    
    log_message "Running GATK BaseRecalibrator for sample: $sample"
    
    # Ensure the known-sites feature file has an accompanying index locally
    ensure_feature_index() {
        local feature_file="$1"
        # Determine expected index based on extension
        if [[ "$feature_file" == *.vcf.gz ]]; then
            if [ ! -f "${feature_file}.tbi" ]; then
                log_message "Index for compressed VCF not found; creating with GATK IndexFeatureFile"
                if ! "${GATK_COMMAND}" IndexFeatureFile --input "${feature_file}"; then
                    error_exit "Failed to create index (.tbi) for ${feature_file}"
                fi
            fi
        elif [[ "$feature_file" == *.vcf ]]; then
            if [ ! -f "${feature_file}.idx" ]; then
                log_message "Index for uncompressed VCF not found; creating with GATK IndexFeatureFile"
                if ! "${GATK_COMMAND}" IndexFeatureFile --input "${feature_file}"; then
                    error_exit "Failed to create index (.idx) for ${feature_file}"
                fi
            fi
        else
            # Other feature types not handled here
            :
        fi
    }
    ensure_feature_index "$known_sites"
    
    if ! "${GATK_COMMAND}" --java-options "-Xmx${memory}" BaseRecalibrator \
        -R "$reference_genome" \
        -I "${sample}_dedup.bam" \
        --known-sites "$known_sites" \
        -O "${sample}_recal.table"; then
        error_exit "GATK BaseRecalibrator failed"
    fi
    
    log_message "GATK BaseRecalibrator completed successfully"
}

# Function to run GATK ApplyBQSR
run_apply_bqsr() {
    local sample="$1"
    local memory="$2"
    local reference_genome="$3"
    
    log_message "Running GATK ApplyBQSR for sample: $sample"
    
    if ! "${GATK_COMMAND}" --java-options "-Xmx${memory}" ApplyBQSR \
        -R "$reference_genome" \
        -I "${sample}_dedup.bam" \
        -bqsr "${sample}_recal.table" \
        -O "${sample}_recal.bam"; then
        error_exit "GATK ApplyBQSR failed"
    fi
    
    log_message "GATK ApplyBQSR completed successfully"
}

# Function to run GATK HaplotypeCaller
run_haplotype_caller() {
    local sample="$1"
    local memory="$2"
    local reference_genome="$3"
    local num_threads="$4"
    
    log_message "Running GATK HaplotypeCaller for sample: $sample"
    
    if ! "${GATK_COMMAND}" --java-options "-Xmx${memory}" HaplotypeCaller \
        -R "$reference_genome" \
        -I "${sample}_recal.bam" \
        -O "${sample}_raw.g.vcf.gz" \
        -ERC GVCF \
        --native-pair-hmm-threads "$num_threads"; then
        error_exit "GATK HaplotypeCaller failed"
    fi
    
    log_message "GATK HaplotypeCaller completed successfully"
}

# Function to run GATK GenotypeGVCFs
run_genotype_gvcfs() {
    local sample="$1"
    local memory="$2"
    local reference_genome="$3"
    
    log_message "Running GATK GenotypeGVCFs for sample: $sample"
    
    if ! "${GATK_COMMAND}" --java-options "-Xmx${memory}" GenotypeGVCFs \
        -R "$reference_genome" \
        -V "${sample}_raw.g.vcf.gz" \
        -O "${sample}_genotyped.vcf.gz"; then
        error_exit "GATK GenotypeGVCFs failed"
    fi
    
    log_message "GATK GenotypeGVCFs completed successfully"
}

# Function: backup current state (sample-specific files only, excludes reference files)
backup_current_state() {
    local step_name="$1"
    local sample="$2"
    local dataset_name="$3"
    
    # Check if backup is enabled
    if [ "${PIPELINE_ENABLE_BACKUP:-true}" != "true" ]; then
        return 0
    fi
    
    local step_num
    case "$step_name" in
        *"Step 2"*|*"Trimmomatic"*)
            step_num=2
            ;;
        *"Step 3"*|*"BWA"*|*"alignment"*)
            step_num=3
            ;;
        *"Step 4"*|*"Duplicate"*|*"MarkDuplicates"*)
            step_num=4
            ;;
        *"Step 5"*|*"BQSR"*|*"recalibration"*)
            step_num=5
            ;;
        *"Step 6"*|*"HaplotypeCaller"*|*"variant calling"*)
            step_num=6
            ;;
        *)
            # Unknown step, skip backup
            return 0
            ;;
    esac
    
    # Allow backups from Step 2 upwards
    if [ "$step_num" -lt 2 ]; then
        return 0
    fi
    
    local backup_dir="${SCRATCH_BASE_PATH%/}/${dataset_name}_backup/${sample}"
    mkdir -p "${backup_dir}"
    
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    log_info "Creating backup for step: $step_name"
    
    if [ "$step_num" -eq 2 ]; then
        # Backup only paired trimmed FASTQs and small logs; skip unpaired and raw FASTQs
        local backed_up_any=false
        if [ -f "${sample}_forward_paired.fastq.gz" ]; then
            rsync -rhivPt "${sample}_forward_paired.fastq.gz" "${backup_dir}/" 2>/dev/null || true
            backed_up_any=true
        fi
        if [ -f "${sample}_reverse_paired.fastq.gz" ]; then
            rsync -rhivPt "${sample}_reverse_paired.fastq.gz" "${backup_dir}/" 2>/dev/null || true
            backed_up_any=true
        fi
        # Optional small logs
        if ls *.log 1>/dev/null 2>&1; then
            rsync -rhivPt *.log "${backup_dir}/" 2>/dev/null || true
            backed_up_any=true
        fi
        if [ "$backed_up_any" = true ]; then
            log_info "Backup completed for step: $step_name (paired trimmed FASTQs and logs)"
        else
            log_info "No paired trimmed FASTQs found to backup for step: $step_name"
        fi
    else
        # Steps 3+ retain existing selective backup (exclude reference files)
        if ls * 2>/dev/null | head -1 >/dev/null 2>&1; then
            rsync -rhivPt \
                --exclude="*.fasta*" \
                --exclude="*.fa" \
                --exclude="*.dict" \
                --exclude="*.fai" \
                --exclude="Known*" \
                --exclude="known*" \
                --exclude="*.tbi" \
                * "${backup_dir}/" 2>/dev/null || true
            log_info "Backup completed for step: $step_name (reference files excluded, sample VCFs included)"
        else
            log_info "No files found to backup for step: $step_name"
        fi
    fi
}

# Function: restore files from backup
restore_from_backup() {
    local sample="$1"
    local dataset_name="$2"
    
    # Check if backup is enabled
    if [ "${PIPELINE_ENABLE_BACKUP:-true}" != "true" ]; then
        return 0
    fi
    
    local backup_dir="${SCRATCH_BASE_PATH%/}/${dataset_name}_backup/${sample}"
    
    if [ -d "$backup_dir" ] && [ "$(ls -A "$backup_dir" 2>/dev/null)" ]; then
        log_info "Restoring files from backup directory: $backup_dir"
        log_info "Files to restore:"
        ls -la "$backup_dir" 2>/dev/null | head -20 || log_info "Could not list backup files"
        
        rsync -rhivPt "$backup_dir"/* . 2>/dev/null || true
        log_info "Restoration completed"
        
        # Verify restoration by listing current directory
        log_info "Files in working directory after restoration:"
        ls -l | head -20
        
        return 0
    else
        log_info "No backup found to restore from"
        return 1
    fi
}

# Function to copy results to RDM
copy_results_to_rdm() {
    local sample="$1"
    local rdm_base_path="$2"
    
    log_message "Copying results to RDM for sample: $sample"
    
    ensure_sample_directories "$rdm_base_path" "$sample"
    
    # Copy BAM files
    if [ -f "${sample}_recal.bam" ]; then
        rsync -rhivPt "${sample}_recal.bam" "${rdm_base_path}/4.BAM/${sample}/" || error_exit "Failed to copy BAM file to RDM"
        rsync -rhivPt "${sample}_recal.bam.bai" "${rdm_base_path}/4.BAM/${sample}/" || error_exit "Failed to copy BAM index to RDM"
    fi
    
    # Copy VCF files
    if [ -f "${sample}_raw.g.vcf.gz" ]; then
        rsync -rhivPt "${sample}_raw.g.vcf.gz" "${rdm_base_path}/5.Individual_VCF/" || error_exit "Failed to copy GVCF file to RDM"
        rsync -rhivPt "${sample}_raw.g.vcf.gz.tbi" "${rdm_base_path}/5.Individual_VCF/" || error_exit "Failed to copy GVCF index to RDM"
    fi
    
    if [ -f "${sample}_genotyped.vcf.gz" ]; then
        rsync -rhivPt "${sample}_genotyped.vcf.gz" "${rdm_base_path}/5.Individual_VCF/" || error_exit "Failed to copy VCF file to RDM"
        rsync -rhivPt "${sample}_genotyped.vcf.gz.tbi" "${rdm_base_path}/5.Individual_VCF/" || error_exit "Failed to copy VCF index to RDM"
    fi
    
    log_message "Results copied to RDM successfully"
}
