#!/bin/bash
# =============================================================================
# STEP 1B FUNCTIONS
# =============================================================================
# Helper functions for Combine GVCFs (Step 1B)
# =============================================================================

# Ensure working directory exists on scratch
ensure_step1b_workdir() {
    local workdir="$1"
    mkdir -p "${workdir}/workspaces"
    mkdir -p "${workdir}/logs"
}

# Utility: setup shared reference genome for array jobs (task 0 only)
setup_shared_reference_genome() {
    local dataset_name="$1"
    local reference_genome="$2"
    
    # Create shared directory
    local shared_base="${SCRATCH_BASE_PATH%/}/${dataset_name}_shared"
    local shared_ref_dir="${shared_base}/Reference_genome"
    local source_dir="${PIPELINE_REFERENCE_DIR:-$(dirname "${reference_genome}")}"
    local ref_basename
    ref_basename="$(basename "${reference_genome}")"
    
    mkdir -p "${shared_ref_dir}"
    
    if [ -f "${shared_ref_dir}/${ref_basename}" ]; then
        log_info "Shared reference genome already exists at ${shared_ref_dir}; verifying assets..."
    else
        log_info "Shared reference genome missing; rsync from ${source_dir}."
        if ! rsync -rhPt "${source_dir%/}/" "${shared_ref_dir}/"; then
            error_exit "Failed to copy reference genome directory to shared location (${source_dir} -> ${shared_ref_dir}). Verify PIPELINE_REFERENCE_DIR."
        fi
    fi
    
    if [ ! -f "${shared_ref_dir}/${ref_basename}" ]; then
        error_exit "Shared reference genome setup incomplete - ${shared_ref_dir}/${ref_basename} missing. Update PIPELINE_REFERENCE_FASTA."
    fi
    
            touch "${shared_base}/.initialized"
    log_info "Shared reference genome is ready at ${shared_ref_dir}"
}

# Utility: get path to shared reference genome if it exists
get_shared_reference_path() {
    local dataset_name="$1"
    local reference_genome="$2"
    
    # Only check for shared files in array job context
    if [ -z "${SLURM_ARRAY_TASK_ID:-}" ]; then
        echo ""
        return 0
    fi
    
    local shared_base="${SCRATCH_BASE_PATH%/}/${dataset_name}_shared"
    local shared_ref_dir="${shared_base}/Reference_genome"
    
    local ref_basename
    ref_basename="$(basename "${reference_genome}")"
    
    # Check if shared reference genome exists
    if [ -f "${shared_ref_dir}/${ref_basename}" ]; then
        echo "${shared_ref_dir}/${ref_basename}"
    fi
}

# Utility: wait for shared reference genome to be initialized (for non-task-0 jobs)
wait_for_shared_reference_genome() {
    local dataset_name="$1"
    local reference_genome="$2"
    
    # Only wait if this is not task 0
    if [ "${SLURM_ARRAY_TASK_ID:-}" = "0" ]; then
        return 0
    fi
    
    local shared_base="${SCRATCH_BASE_PATH%/}/${dataset_name}_shared"
    local marker_file="${shared_base}/.initialized"
    local timeout="${PIPELINE_SHARED_REF_TIMEOUT:-120}"
    local poll_interval="${PIPELINE_SHARED_REF_POLL_INTERVAL:-3}"
    local elapsed=0
    
    log_info "Waiting for shared reference genome initialization (timeout: ${timeout}s, polling every ${poll_interval}s)"
    
    while [ $elapsed -lt $timeout ]; do
        # Check marker file exists
        if [ -f "${marker_file}" ]; then
            # Verify actual reference genome exists (not just marker)
            local ref_basename
            ref_basename="$(basename "${reference_genome}")"
            
            if [ -f "${shared_base}/Reference_genome/${ref_basename}" ]; then
                log_info "Shared reference genome is ready (waited ${elapsed}s)"
                return 0
            fi
        fi
        
        sleep "$poll_interval"
        elapsed=$((elapsed + poll_interval))
        
        # Log progress every 15 seconds
        if [ $((elapsed % 15)) -eq 0 ]; then
            log_info "Still waiting for shared reference genome... (${elapsed}/${timeout}s)"
        fi
    done
    
    error_exit "Timeout waiting for shared reference genome initialization (${timeout}s exceeded)"
}

# Copy reference genome files to working directory (matching original script approach)
# Uses symlinks if shared reference genome exists, otherwise copies directly
copy_reference_genome_to_workdir() {
    local reference_genome="$1"
    local workdir="$2"
    local dataset_name="${3:-}"
    
    # Check if shared reference genome exists (for array jobs)
    local shared_ref_path
    if [ -n "${dataset_name}" ] && [ -n "${SLURM_ARRAY_TASK_ID:-}" ]; then
        shared_ref_path="$(get_shared_reference_path "${dataset_name}" "${reference_genome}")"
    else
        shared_ref_path=""
    fi
    
    local ref_basename
    ref_basename="$(basename "${reference_genome}")"
    local local_ref="${workdir}/${ref_basename}"
    
    if [ -n "${shared_ref_path}" ]; then
        # Use symlinks for shared reference genome
        log_info "Using shared reference genome (creating symlinks)"
        
        # Create symlinks
        ln -sf "${shared_ref_path}" "${local_ref}" || error_exit "Failed to create symlink for reference genome"
        
        # Create symlinks for index files if they exist
        if [ -f "${shared_ref_path}.fai" ]; then
            ln -sf "${shared_ref_path}.fai" "${local_ref}.fai" || log_info "Warning: Failed to create symlink for reference genome index (.fai), continuing anyway"
        fi
        
        local ref_dict="${shared_ref_path%.*}.dict"
        if [ -f "${ref_dict}" ]; then
            local ref_dict_basename
            ref_dict_basename="$(basename "${ref_dict}")"
            ln -sf "${ref_dict}" "${workdir}/${ref_dict_basename}" || log_info "Warning: Failed to create symlink for reference genome dictionary (.dict), continuing anyway"
        fi
        
        log_info "Reference genome files linked successfully"
    else
        # Copy files directly (non-array job or shared files not available)
        log_info "Copying reference genome files to working directory"
        
        rsync -rhivPt "${reference_genome}" "${local_ref}" || error_exit "Failed to copy reference genome to working directory"
        
        # Copy reference index (.fai) if it exists
        if [ -f "${reference_genome}.fai" ]; then
            rsync -rhivPt "${reference_genome}.fai" "${local_ref}.fai" || log_info "Warning: Failed to copy reference genome index (.fai), continuing anyway"
        fi
        
        # Copy reference dictionary (.dict) if it exists
        local ref_dict="${reference_genome%.*}.dict"
        if [ -f "${ref_dict}" ]; then
            local ref_dict_basename
            ref_dict_basename="$(basename "${ref_dict}")"
            rsync -rhivPt "${ref_dict}" "${workdir}/${ref_dict_basename}" || log_info "Warning: Failed to copy reference genome dictionary (.dict), continuing anyway"
        fi
        
        log_info "Reference genome files copied successfully"
    fi
    
    # Return local reference path
    echo "${local_ref}"
}

# Build sample map referencing GVCF files on RDM storage
build_sample_map() {
    local sample_map_path="$1"
    local rdm_base_path="$2"

    log_info "Creating sample map at ${sample_map_path}"
    > "${sample_map_path}"

    local gvcf_dir="${rdm_base_path}/5.Individual_VCF"
    local gvcfs=("${gvcf_dir}"/*_raw.g.vcf.gz)

    if [ ${#gvcfs[@]} -eq 0 ]; then
        error_exit "No *_raw.g.vcf.gz files found in ${gvcf_dir}"
    fi

    for gvcf_file in "${gvcfs[@]}"; do
        [ -f "${gvcf_file}" ] || continue
        local sample_name
        sample_name=$(basename "${gvcf_file}" "_raw.g.vcf.gz")
        printf "%s\t%s\n" "${sample_name}" "${gvcf_file}" >> "${sample_map_path}"
    done

    if [ ! -s "${sample_map_path}" ]; then
        error_exit "Sample map ${sample_map_path} is empty"
    fi

    log_info "Sample map created with $(wc -l < "${sample_map_path}") entries"
}

# Run GenomicsDBImport for a chromosome/interval
run_genomics_db_import() {
    local chromosome="$1"
    local workdir="$2"
    local reference_genome="$3"
    local sample_map="$4"
    local memory="$5"

    log_info "Running GenomicsDBImport for ${chromosome}"

    local db_workspace="${workdir}/workspaces/genomicsdb_${chromosome}"
    if [ -d "${db_workspace}" ]; then
        log_warn "Removing existing GenomicsDB workspace for ${chromosome}: ${db_workspace}"
        rm -rf "${db_workspace}"
    fi
    mkdir -p "${workdir}/tmp"

    if ! "${GATK_COMMAND}" --java-options "-Xmx${memory}" GenomicsDBImport \
        --genomicsdb-workspace-path "${db_workspace}" \
        --batch-size "${GATK_BATCH_SIZE}" \
        --sample-name-map "${sample_map}" \
        --intervals "${chromosome}" \
        --reference "${reference_genome}" \
        --reader-threads "${GATK_READER_THREADS}" \
        --genomicsdb-shared-posixfs-optimizations true \
        --merge-input-intervals true \
        --consolidate true \
        --overwrite-existing-genomicsdb-workspace true \
        --tmp-dir "${workdir}/tmp"; then
        error_exit "GenomicsDBImport failed for ${chromosome}"
    fi

    log_info "GenomicsDBImport completed for ${chromosome}"
}

# Run GenotypeGVCFs using a GenomicsDB workspace
run_genotype_gvcfs() {
    local chromosome="$1"
    local workdir="$2"
    local reference_genome="$3"
    local memory="$4"
    local output_file="$5"

    log_info "Running GenotypeGVCFs for ${chromosome}"

    if ! "${GATK_COMMAND}" --java-options "-Xmx${memory}" GenotypeGVCFs \
        -R "${reference_genome}" \
        -V "gendb://${workdir}/workspaces/genomicsdb_${chromosome}" \
        -O "${output_file}"; then
        error_exit "GenotypeGVCFs failed for ${chromosome}"
    fi

    log_info "GenotypeGVCFs completed for ${chromosome}"
}

# Copy consolidated VCF (and index if available) back to RDM
copy_consolidated_vcf() {
    local output_file="$1"
    local rdm_base_path="$2"
    local chromosome="$3"

    local destination="${rdm_base_path}/7.Consolidated_VCF"
    mkdir -p "${destination}"

    log_info "Copying ${chromosome} consolidated VCF to ${destination}"

    rsync -rhivPt "${output_file}" "${destination}/" || error_exit "Failed to copy VCF for ${chromosome}"
    if [ -f "${output_file}.tbi" ]; then
        rsync -rhivPt "${output_file}.tbi" "${destination}/" || error_exit "Failed to copy VCF index for ${chromosome}"
    fi
}

# Remove GenomicsDB workspace for a chromosome to free scratch space
cleanup_chromosome_workspace() {
    local workdir="$1"
    local chromosome="$2"
    rm -rf "${workdir}/workspaces/genomicsdb_${chromosome}"
}

# Obtain chromosome list from reference
get_chromosome_list() {
    local reference_genome="$1"

    if [ "${#PIPELINE_CHROMOSOMES[@]}" -gt 0 ]; then
        log_info "Using configured chromosome list (${#PIPELINE_CHROMOSOMES[@]} entries) for ${reference_genome}"
        printf '%s\n' "${PIPELINE_CHROMOSOMES[@]}"
        return 0
    fi

    local fai_path="${reference_genome}.fai"
    if [ -f "${fai_path}" ]; then
        log_warn "PIPELINE_CHROMOSOME_LIST empty; deriving contigs from ${fai_path}"
        awk '{print $1}' "${fai_path}"
        return 0
    fi

    local reference_dir="${PIPELINE_REFERENCE_DIR:-$(dirname "${reference_genome}")}"
    fai_path="${reference_dir%/}/$(basename "${reference_genome}").fai"
    if [ -f "${fai_path}" ]; then
        log_warn "Using FASTA index from reference directory: ${fai_path}"
        awk '{print $1}' "${fai_path}"
        return 0
    fi

    error_exit "No chromosomes configured and unable to locate $(basename "${reference_genome}").fai. Copy the Reference_Genome assets (see Data_management) or set PIPELINE_CHROMOSOME_LIST."
}
