#!/bin/bash
# =============================================================================
# PROFESSIONAL LOGGING SYSTEM FOR GATK PIPELINE
# =============================================================================
# This logging system works within SLURM's logging framework
# It provides structured logging for all pipeline components

# Logging configuration
LOG_LEVELS=("DEBUG" "INFO" "WARN" "ERROR" "FATAL")
DEFAULT_LOG_LEVEL="INFO"
LOG_TO_CONSOLE=true
LOG_TO_FILE=true

# Colors for better user experience
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
PURPLE='\033[0;35m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

_default_scratch_base() {
    if [ -n "${SCRATCH_BASE_PATH:-}" ]; then
        echo "${SCRATCH_BASE_PATH%/}"
        return
    fi
    if [ -d "/scratch/user/${USER:-}" ]; then
        echo "/scratch/user/${USER}"
        return
    fi
    if [ -d "/scratch/${USER:-}" ]; then
        echo "/scratch/${USER}"
        return
    fi
    if [ -d "/scratch" ]; then
        echo "/scratch"
        return
    fi
    echo "/tmp"
}

resolve_log_root() {
    local dataset_name="$1"
    local subdir="${2:-pipeline}"
    local base_root=""

    if [ -n "${PIPELINE_LOG_DIR_OVERRIDE:-}" ]; then
        base_root="${PIPELINE_LOG_DIR_OVERRIDE%/}"
    elif [ -n "${LOG_BASE_PATH:-}" ]; then
        base_root="${LOG_BASE_PATH%/}"
    else
        local scratch_base
        scratch_base="$(_default_scratch_base)"
        base_root="${scratch_base%/}/logs"
        LOG_BASE_PATH="${base_root}"
        export LOG_BASE_PATH
    fi

    if [[ "${base_root}" != /* ]]; then
        local current_dir
        current_dir="$(pwd 2>/dev/null || echo '/tmp')"
        base_root="${current_dir%/}/${base_root#./}"
    fi

    local root="${base_root}"
    if [ -n "${dataset_name}" ] && [ -z "${PIPELINE_LOG_DIR_OVERRIDE:-}" ]; then
        root="${root%/}/${dataset_name}"
    fi

    if [ -n "${subdir}" ]; then
        root="${root%/}/${subdir}"
    fi

    if ! mkdir -p "${root}" 2>/dev/null; then
        echo "Warning: Failed to create log directory: ${root}, using /tmp" >&2
        root="/tmp"
        mkdir -p "${root}" 2>/dev/null || true
    else
        root="$(cd "${root}" 2>/dev/null && pwd)" || root="${root}"
    fi

    printf '%s\n' "${root}"
}

# Initialize logging for any module
init_logging() {
    local module_name="$1"
    local log_type="$2"  # "wrapper", "pipeline", "slurm", "job"
    local dataset_name="${3:-${DATASET_NAME:-${DATASET:-}}}"
    local log_subdir="${4:-}"
    local timestamp
    timestamp=$(date '+%Y%m%d_%H%M%S')

    case "$log_type" in
        "pipeline"|"job"|"wrapper")
            log_subdir="${log_subdir:-pipeline}"
            ;;
        "slurm")
            log_subdir="${log_subdir:-slurm}"
            ;;
        *)
            log_subdir="${log_subdir:-pipeline}"
            ;;
    esac

    LOG_DIR="$(resolve_log_root "${dataset_name}" "${log_subdir}")"

    case "$log_type" in
        "wrapper")
            LOG_FILE="${module_name}_wrapper_${timestamp}.log"
            ;;
        "pipeline"|"job")
            LOG_FILE="pipeline_${module_name}_${timestamp}.log"
            ;;
        "slurm")
            LOG_FILE="slurm_${module_name}_${timestamp}.log"
            ;;
        *)
            LOG_FILE="${module_name}_${timestamp}.log"
            ;;
    esac

    LOG_FILE="${LOG_DIR%/}/${LOG_FILE}"

    LOG_LEVEL="${LOG_LEVEL:-$DEFAULT_LOG_LEVEL}"
    MODULE_NAME="$module_name"

    log_info "Logging system initialized for module: $module_name"
    log_info "Log type: $log_type"
    log_info "Log directory: $LOG_DIR"
    log_info "Log file: $LOG_FILE"
    log_info "Log level: $LOG_LEVEL"
}

# Get numeric log level
get_log_level_num() {
    local level="$1"
    case "$level" in
        "DEBUG") echo 0;;
        "INFO")  echo 1;;
        "WARN")  echo 2;;
        "ERROR") echo 3;;
        "FATAL") echo 4;;
        *)       echo 1;; # Default to INFO
    esac
}

# Core logging function
log_message() {
    local level="${1:-INFO}"
    local message="${2:-}"
    local module="${MODULE_NAME:-unknown}"
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    local log_entry="[$timestamp] [$level] [$module] $message"
    
    # Check if we should log this level
    local effective_log_level="${LOG_LEVEL:-$DEFAULT_LOG_LEVEL}"
    local current_level_num=$(get_log_level_num "$effective_log_level")
    local message_level_num=$(get_log_level_num "$level")
    
    # Use absolute path for log file, fallback to /tmp if not set
    local log_file="${LOG_FILE:-/tmp/gatk_pipeline.sh.log}"
    
    # Ensure log_file is absolute
    if [[ "$log_file" != /* ]]; then
        local current_dir
        current_dir="$(pwd 2>/dev/null || echo '/tmp')"
        log_file="${current_dir}/${log_file}"
    fi
    
    if [ $message_level_num -ge $current_level_num ]; then
        # Log to file (ensure directory exists)
        if [ "$LOG_TO_FILE" = "true" ]; then
            local log_file_dir
            log_file_dir="$(dirname "$log_file")"
            mkdir -p "$log_file_dir" 2>/dev/null || true
            echo "$log_entry" >> "$log_file" 2>/dev/null || {
                # If writing fails, try /tmp as fallback
                echo "$log_entry" >> "/tmp/gatk_pipeline.sh_${MODULE_NAME:-unknown}.log" 2>/dev/null || true
            }
        fi
        
        # Log to console with colors
        if [ "$LOG_TO_CONSOLE" = "true" ]; then
            case "$level" in
                "DEBUG") >&2 echo -e "${PURPLE}$log_entry${NC}";;
                "INFO")  >&2 echo -e "${BLUE}$log_entry${NC}";;
                "WARN")  >&2 echo -e "${YELLOW}$log_entry${NC}";;
                "ERROR") >&2 echo -e "${RED}$log_entry${NC}";;
                "FATAL") >&2 echo -e "${RED}$log_entry${NC}";;
                *)       >&2 echo "$log_entry";;
            esac
        fi
    fi
}

# Convenience logging functions
log_debug() { log_message "DEBUG" "$1"; }
log_info()  { log_message "INFO"  "$1"; }
log_warn()  { log_message "WARN"  "$1"; }
log_error() { log_message "ERROR" "$1"; }
log_fatal() { log_message "FATAL" "$1"; }

# Module-specific logging functions
log_module_info() { log_message "INFO" "[$MODULE_NAME] $1"; }
log_module_error() { log_message "ERROR" "[$MODULE_NAME] $1"; }
log_module_warn() { log_message "WARN" "[$MODULE_NAME] $1"; }
log_module_debug() { log_message "DEBUG" "[$MODULE_NAME] $1"; }

# SLURM-aware logging functions
log_slurm_submission() {
    local job_id="$1"
    local slurm_script="$2"
    local dataset_name="$3"
    local sample_count="$4"
    local slurm_log_dir="${LOG_BASE_PATH}/${dataset_name}"
    mkdir -p "${slurm_log_dir}"
    
    log_info "SLURM job submitted successfully"
    log_info "Job ID: $job_id"
    log_info "SLURM script: $slurm_script"
    log_info "Sample count: $sample_count"
    log_info "SLURM log directory: $slurm_log_dir"
    log_info "Monitor with: squeue -u \$USER"
    log_info "Cancel with: scancel $job_id"
    log_info "View SLURM logs: ls $slurm_log_dir/"
    
    # Create a job tracking file for easy reference
    local job_tracker="${slurm_log_dir}/job_${job_id}.info"
    cat > "$job_tracker" << EOF
# SLURM Job Information
Job ID: $job_id
Dataset: $dataset_name
Sample Count: $sample_count
SLURM Script: $slurm_script
Submitted: $(date '+%Y-%m-%d %H:%M:%S')
SLURM Log Directory: $slurm_log_dir
Wrapper Log: $LOG_FILE

# Monitoring Commands
Monitor: squeue -u \$USER
Cancel: scancel $job_id
View Logs: ls $slurm_log_dir/
Job Details: scontrol show job $job_id
EOF
    
    log_info "Job tracking file created: $job_tracker"
}

log_slurm_status() {
    local job_id="$1"
    local status="$2"
    
    case "$status" in
        "RUNNING")
            log_info "SLURM job $job_id is running"
            ;;
        "COMPLETED")
            log_info "SLURM job $job_id completed successfully"
            ;;
        "FAILED")
            log_error "SLURM job $job_id failed"
            ;;
        "CANCELLED")
            log_warn "SLURM job $job_id was cancelled"
            ;;
        *)
            log_info "SLURM job $job_id status: $status"
            ;;
    esac
}

# Function to print colored messages (legacy support)
print_message() {
    local color=$1
    local message=$2
    echo -e "${color}${message}${NC}"
    log_info "$message"  # Also log the message
}

# Error exit function - logs error and exits with code 1
error_exit() {
    local error_message="$1"
    log_fatal "$error_message"
    exit 1
}
