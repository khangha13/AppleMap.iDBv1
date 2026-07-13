#!/bin/bash
# =============================================================================
# SLURM JOB MANAGEMENT UTILITIES FOR GATK PIPELINE
# =============================================================================
# SLURM job management functions used across all pipeline steps

build_slurm_export_arg() {
    local export_vars=(
        PIPELINE_ROOT
        PIPELINE_ENV_FILE
        PIPELINE_DIR_NAME
        PIPELINE_HOME_CANDIDATE
        PIPELINE_LOG_DIR_OVERRIDE
        PIPELINE_RDM_BASE
        PIPELINE_RDM_DATASETS
        PIPELINE_SCRATCH_BASE
        PIPELINE_LOG_BASE
        PIPELINE_REFERENCE_DIR
        PIPELINE_REFERENCE_FASTA
        PIPELINE_KNOWN_SITES_VCF
        PIPELINE_ADAPTER_FASTA
        PIPELINE_SLURM_ACCOUNT
        PIPELINE_SLURM_PARTITION
        PIPELINE_SLURM_QOS
        PIPELINE_SLURM_NODES
        PIPELINE_SLURM_NTASKS
        PIPELINE_SLURM_SCRIPT_DIR
        PIPELINE_WORK_DIR
        PIPELINE_LOG_LEVEL
        PIPELINE_VERBOSE_LOGGING
        PIPELINE_GLOBAL_DRY_RUN
        PIPELINE_ENABLE_BACKUP
        PIPELINE_BACKUP_STEPS
        PIPELINE_SHARED_REF_TIMEOUT
        PIPELINE_SHARED_REF_POLL_INTERVAL
        PIPELINE_CHROMOSOME_LIST
        RDM_BASE_PATH
        RDM_DATASETS_PATH
        SCRATCH_BASE_PATH
        LOG_BASE_PATH
        MASTER_LOG_DIR
        MASTER_LOG_FILE
        WRAPPER_LOG_PATH
        DEFAULT_LOG_LEVEL
        LOG_TO_CONSOLE
        LOG_TO_FILE
        DEBUG_MODE
        TEST_MODE
        REFERENCE_GENOME
        REFERENCE_GENOME_INDEX
        REFERENCE_GENOME_DICT
        REFERENCE_GENOME_ORIGINAL
        KNOWN_SITES_VCF
        KNOWN_SITES_INDEX
        ADAPTER_FILE
        STEP1A_CPUS
        STEP1A_MEMORY
        STEP1A_TIME
        STEP1A_ARRAY_LIMIT
        STEP1A_SAMTOOLS_MODULE
        STEP1B_CPUS
        STEP1B_MEMORY
        STEP1B_TIME
        STEP1B_ARRAY_LIMIT
        STEP1B_CHROMOSOME_LIST
        STEP1B_FAILURE_FLAG_PATH
        STEP1B_RUNNING_FLAG_PATH
        STEP1B_FAILURE_CONTEXT
        STEP1B_GATK_MODULE
        STEP1C_CPUS
        STEP1C_MEMORY
        STEP1C_TIME
        STEP1C_ACCOUNT
        STEP1C_PARTITION
        STEP1C_NODES
        STEP1C_NTASKS
        STEP1C_SELF_IMPUTE
        STEP1C_ARRAY_LIMIT
        STEP1C_CONSTRAINT
        STEP1C_DEBUG_BEAGLE
        STEP1C_PERSIST_ARTIFACTS
        STEP1D_CPUS
        STEP1D_MEMORY
        STEP1D_TIME
        STEP1D_ACCOUNT
        STEP1D_PARTITION
        STEP1D_NODES
        STEP1D_NTASKS
        STEP1D_CONSTRAINT
        STEP1D_PCA_DIR
        STEP1D_PCA_MERGED_PATTERN
        STEP1D_PCA_FORCE_CONCAT
        STEP1D_PCA_MERGED_EXCLUDE_CHR
        STEP1D_CACHE_DIR
        STEP1D_PACKAGE_DIR
        STEP1D_PARQUET_COMPRESSION
        GATK_COMMAND
        GATK_BATCH_SIZE
        GATK_READER_THREADS
        BWA_MEM_FLAGS
        BWA_PATH
        SAMTOOLS_PATH
        SAMTOOLS_MODULE
        FASTQC_PATH
        TRIMMOMATIC_JAR
        BCFTOOLS_BIN
        PLINK2_BIN
        BEAGLE_JAR
        PYTHON_BIN
    )
    local export_items=()
    local var value
    for var in "${export_vars[@]}"; do
        value="${!var:-}"
        # Slurm --export is comma-delimited; comma-valued config is reloaded
        # from pipeline_config.sh in the child job instead of being inlined here.
        if [[ "${value}" == *,* ]]; then
            continue
        fi
        if [ -n "${value}" ]; then
            export_items+=("${var}=${value}")
        fi
    done

    if [ "${#export_items[@]}" -eq 0 ]; then
        echo "--export=NONE"
    else
        local IFS=,
        echo "--export=${export_items[*]}"
    fi
}

# Function to create SLURM script from template
create_slurm_script() {
    local template="$1"
    local config="$2"
    local output="$3"
    local dataset_name="$4"
    local step_name="$5"
    local pipeline_root="$6"
    
    log_info "Creating SLURM script: $output"

    if [ ! -f "$template" ]; then
        log_error "SLURM template not found: $template"
        return 1
    fi

    if [ -z "$pipeline_root" ]; then
        log_error "PIPELINE_ROOT not provided to create_slurm_script"
        return 1
    fi

    local output_dir
    output_dir="$(dirname "$output")"
    if ! mkdir -p "$output_dir"; then
        log_error "Failed to create output directory: $output_dir"
        return 1
    fi
    
    if [ ! -d "$output_dir" ]; then
        log_error "Output directory does not exist after creation: $output_dir"
        return 1
    fi
    
    # Parse configuration using safer method
    declare -A config_map=()
    while IFS='=' read -r key value; do
        [ -z "$key" ] && continue
        # Handle values that may contain '=' by taking everything after first '='
        if [[ "$key" == *"="* ]]; then
            # This shouldn't happen with proper config format, but handle it
            continue
        fi
        config_map["$key"]="$value"
    done <<< "$config"
    
    # Extract configuration values with defaults
    local job_name="${config_map[job_name]:-}"
    local account="${config_map[account]:-}"
    local partition="${config_map[partition]:-}"
    local nodes="${config_map[nodes]:-}"
    local ntasks="${config_map[ntasks]:-}"
    local cpus_per_task="${config_map[cpus_per_task]:-}"
    local time_limit="${config_map[time_limit]:-}"
    local memory="${config_map[memory]:-}"
    local array_max="${config_map[array_max]:-}"
    local qos="${config_map[qos]:-}"
    local constraint="${config_map[constraint]:-}"
    
    # Validate required fields
    local missing_fields=()
    [ -z "$job_name" ] && missing_fields+=("job_name")
    [ -z "$account" ] && missing_fields+=("account")
    [ -z "$partition" ] && missing_fields+=("partition")
    [ -z "$nodes" ] && missing_fields+=("nodes")
    [ -z "$ntasks" ] && missing_fields+=("ntasks")
    [ -z "$cpus_per_task" ] && missing_fields+=("cpus_per_task")
    [ -z "$time_limit" ] && missing_fields+=("time_limit")
    [ -z "$memory" ] && missing_fields+=("memory")
    [ -z "$array_max" ] && missing_fields+=("array_max")
    
    if [ ${#missing_fields[@]} -gt 0 ]; then
        log_error "Missing required SLURM config fields: ${missing_fields[*]}"
        return 1
    fi
    
    local log_root=""
    if command -v resolve_log_root >/dev/null 2>&1; then
        log_root="$(resolve_log_root "${dataset_name}" "slurm")"
    else
        log_root="${LOG_BASE_PATH%/}/${dataset_name}/slurm"
        mkdir -p "${log_root}" 2>/dev/null || log_root="/tmp"
    fi

    # Create SLURM script header
    if ! cat > "$output" << EOF
#!/bin/bash -l
#SBATCH --job-name=${job_name}
#SBATCH --ntasks=${ntasks}
#SBATCH --account=${account}
#SBATCH --partition=${partition}
#SBATCH --nodes=${nodes}
#SBATCH --cpus-per-task=${cpus_per_task}
#SBATCH --time=${time_limit}
#SBATCH --mem=${memory}
#SBATCH --array=0-${array_max}
#SBATCH --chdir=${log_root}
#SBATCH -o ${log_root}/${step_name}_%A_%a_%x_%j.output
#SBATCH -e ${log_root}/${step_name}_%A_%a_%x_%j.error

# =============================================================================
# SLURM-GENERATED SCRIPT FOR APPLE GATK PIPELINE - ${step_name^^}
# =============================================================================
# Generated by bin/gatk_pipeline.sh on $(date)
# Dataset: ${dataset_name}
# Array Range: 0-${array_max}
# SLURM Config: ${cpus_per_task} CPUs, ${memory} RAM, ${time_limit} time
# =============================================================================

EOF
    then
        log_error "Failed to write SLURM script header to: $output"
        return 1
    fi

    # Append constraint line if provided
    if [ -n "${constraint}" ]; then
        echo "#SBATCH --constraint=${constraint}" >> "$output"
        echo "" >> "$output"
    fi
    
    local py_bin="${PYTHON_BIN:-python3}"
    if ! command -v "$py_bin" >/dev/null 2>&1; then
        py_bin="python"
    fi
    # Use stdin ("-") to avoid Python treating the template path as a script file
    if ! "$py_bin" - <<'PY' "$template" "$output" "$pipeline_root"; then
import sys
template, output, pipeline_root = sys.argv[1:4]
with open(template, encoding="utf-8") as src:
    lines = src.readlines()
with open(output, "a", encoding="utf-8") as dst:
    for line in lines[1:]:
        if line.startswith("PIPELINE_ROOT=") and "BASH_SOURCE" in line:
            dst.write(f'PIPELINE_ROOT="{pipeline_root}"\n')
        else:
            dst.write(line.replace("__PIPELINE_ROOT_PLACEHOLDER__", pipeline_root))
PY
        log_error "Failed to append template body from $template to $output"
        return 1
    fi
    
    # Verify the file was created successfully
    if [ ! -f "$output" ]; then
        log_error "SLURM script file was not created: $output"
        return 1
    fi
    
    if [ ! -s "$output" ]; then
        log_error "SLURM script file is empty: $output"
        return 1
    fi
    
    # Make it executable
    if ! chmod +x "$output"; then
        log_error "Failed to make script executable: $output"
        return 1
    fi
    
    log_info "✓ Created SLURM script: $output"
    log_info "  Job name: $job_name"
    log_info "  Array range: 0-${array_max}"
    log_info "  SLURM config: ${cpus_per_task} CPUs, ${memory} RAM, ${time_limit} time"
}

# Function to submit job with logging
submit_job() {
    local script="$1"
    local params="$2"
    local dataset_name="$3"
    local step_name="$4"
    
    log_info "Submitting $step_name job: $script with parameters: $params"

    # Note: File existence check removed - create_slurm_script() already verifies
    # the file exists before returning. On network filesystems, immediate checks
    # can fail due to metadata caching even though the file exists.

    # Create log directory
    local submit_log_root=""
    if command -v resolve_log_root >/dev/null 2>&1; then
        submit_log_root="$(resolve_log_root "${dataset_name}" "slurm")"
    else
        submit_log_root="${LOG_BASE_PATH%/}/${dataset_name}/slurm"
        mkdir -p "${submit_log_root}" 2>/dev/null || submit_log_root="/tmp"
    fi
    log_info "SLURM log directory: ${submit_log_root}"
    
    # Submit the job
    # Resolve script to absolute path to avoid path resolution issues
    # Strip any trailing newlines and whitespace that might come from command substitution
    local script_clean
    script_clean=$(printf '%s' "$script" | tr -d '\n\r' | sed 's/[[:space:]]*$//' | head -n1)
    
    # Get absolute path - check if already absolute first
    local script_abs
    if [[ "$script_clean" = /* ]]; then
        # Already absolute path - use it directly
        script_abs="$script_clean"
    else
        # Relative path - resolve it
        local script_dir script_file
        script_dir="$(cd "$(dirname "$script_clean")" 2>/dev/null && pwd)"
        if [ $? -ne 0 ] || [ -z "$script_dir" ]; then
            log_error "Failed to resolve script directory for: $script_clean"
            return 1
        fi
        script_file="$(basename "$script_clean")"
        script_abs="${script_dir}/${script_file}"
    fi
    
    local export_arg
    export_arg="$(build_slurm_export_arg)"
    
    local sbatch_output
    if ! sbatch_output=$(sbatch --parsable "${export_arg}" "$script_abs" $params 2>&1); then
        log_error "Failed to submit $step_name job via sbatch: ${sbatch_output}"
        log_error "Script path attempted: ${script_abs}"
        return 1
    fi
    local job_id
    job_id="$(printf '%s\n' "$sbatch_output" | head -n1 | cut -d';' -f1 | tr -d '\r')"
    if [ -z "${job_id}" ]; then
        log_error "sbatch returned empty job ID (raw output: ${sbatch_output})"
        return 1
    fi
    
    log_info "$step_name job submitted successfully - Job ID: $job_id"
    echo "$job_id"
}

# Function to monitor job status
monitor_job() {
    local job_id="$1"
    local timeout="${2:-3600}"  # Default 1 hour timeout
    
    log_info "Monitoring job: $job_id"
    
    local start_time=$(date +%s)
    local end_time=$((start_time + timeout))
    
    while [ $(date +%s) -lt $end_time ]; do
        if ! squeue -j "$job_id" >/dev/null 2>&1; then
            log_info "Job $job_id completed"
            return 0
        fi
        
        log_debug "Job $job_id still running..."
        sleep 60  # Check every minute
    done
    
    log_warn "Job $job_id monitoring timed out after $timeout seconds"
    return 1
}

# Function to get job status
get_job_status() {
    local job_id="$1"
    
    if ! squeue -j "$job_id" >/dev/null 2>&1; then
        echo "COMPLETED"
    else
        local status=$(squeue -j "$job_id" --noheader --format="%.2T" 2>/dev/null)
        echo "${status:-UNKNOWN}"
    fi
}

# Function to cancel job
cancel_job() {
    local job_id="$1"
    
    log_info "Cancelling job: $job_id"
    
    if scancel "$job_id"; then
        log_info "Job $job_id cancelled successfully"
        return 0
    else
        log_error "Failed to cancel job $job_id"
        return 1
    fi
}

# Function to get job information
get_job_info() {
    local job_id="$1"
    
    log_info "Getting job information for: $job_id"
    
    if scontrol show job "$job_id" >/dev/null 2>&1; then
        scontrol show job "$job_id"
    else
        log_warn "Job $job_id not found or completed"
        return 1
    fi
}

# Function to get user input with validation
get_user_input() {
    local prompt=$1
    local default=$2
    local input
    
    if [ -n "$default" ]; then
        read -p "$prompt [$default]: " input
        input=${input:-$default}
    else
        read -p "$prompt: " input
    fi
    
    echo "$input"
}

# Function to confirm before proceeding
confirm_action() {
    local message=$1
    local response
    
    while true; do
        read -p "$message (y/n): " response
        case $response in
            [Yy]* ) return 0;;
            [Nn]* ) return 1;;
            * ) echo "Please answer yes or no.";;
        esac
    done
}
