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

# Initialize logging for any module
init_logging() {
    local module_name="$1"
    local log_type="$2"  # "wrapper", "pipeline", "slurm"
    local timestamp=$(date '+%Y%m%d_%H%M%S')
    
    # Determine log directory using absolute paths
    local log_dir_base=""
    
    # Allow explicit override for special cases (e.g., master script run folders)
    if [ -n "${PIPELINE_LOG_DIR_OVERRIDE:-}" ]; then
        log_dir_base="${PIPELINE_LOG_DIR_OVERRIDE}"
    # For pipeline logs in SLURM context, use LOG_BASE_PATH if available
    elif [ "$log_type" = "pipeline" ] && [ -n "${LOG_BASE_PATH:-}" ]; then
        # Extract dataset name from DATASET_PATH if available (for SLURM array jobs)
        if [ -n "${DATASET_PATH:-}" ] && [ -d "${DATASET_PATH}" ]; then
            local dataset_name=$(basename "${DATASET_PATH}")
            log_dir_base="${LOG_BASE_PATH}/${dataset_name}"
        else
            log_dir_base="${LOG_BASE_PATH}"
        fi
    # For wrapper logs, use WRAPPER_LOG_PATH if available, otherwise relative to current dir
    elif [ "$log_type" = "wrapper" ] && [ -n "${WRAPPER_LOG_PATH:-}" ]; then
        log_dir_base="${WRAPPER_LOG_PATH}"
    # For other cases, use absolute path based on current working directory
    else
        local current_dir
        current_dir="$(pwd 2>/dev/null || echo '/tmp')"
        case "$log_type" in
            "wrapper")
                log_dir_base="${current_dir}/logs/wrapper"
                ;;
            "pipeline")
                log_dir_base="${current_dir}/logs/pipeline"
                ;;
            "slurm")
                log_dir_base="${current_dir}/logs/slurm"
                ;;
            *)
                log_dir_base="${current_dir}/logs"
                ;;
        esac
    fi
    
    # Set log file name
    case "$log_type" in
        "wrapper")
            LOG_FILE="${module_name}_wrapper_${timestamp}.log"
            ;;
        "pipeline")
            LOG_FILE="pipeline_${module_name}_${timestamp}.log"
            ;;
        "slurm")
            LOG_FILE="slurm_${module_name}_${timestamp}.log"
            ;;
        *)
            LOG_FILE="${module_name}_${timestamp}.log"
            ;;
    esac
    
    # Ensure log directory is absolute
    if [[ "$log_dir_base" != /* ]]; then
        # Convert relative path to absolute
        local current_dir
        current_dir="$(pwd 2>/dev/null || echo '/tmp')"
        # Remove leading ./ if present
        log_dir_base="${log_dir_base#./}"
        LOG_DIR="${current_dir}/${log_dir_base}"
    else
        LOG_DIR="${log_dir_base}"
    fi
    
    # Normalize the path by resolving parent directory if it exists
    local parent_dir
    parent_dir="$(dirname "${LOG_DIR}")"
    if [ -d "$parent_dir" ]; then
        # Parent exists, normalize it
        parent_dir="$(cd "$parent_dir" 2>/dev/null && pwd)" || parent_dir="$(dirname "${LOG_DIR}")"
        LOG_DIR="${parent_dir}/$(basename "${LOG_DIR}")"
    fi
    # If parent doesn't exist, we'll create it with mkdir -p below
    
    # Create log directory if it doesn't exist (with error handling)
    if ! mkdir -p "$LOG_DIR" 2>/dev/null; then
        # Fallback to /tmp if directory creation fails
        echo "Warning: Failed to create log directory: $LOG_DIR, using /tmp" >&2
        LOG_DIR="/tmp"
        mkdir -p "$LOG_DIR" 2>/dev/null || true
    else
        # After creation, normalize the path
        LOG_DIR="$(cd "$LOG_DIR" 2>/dev/null && pwd)" || LOG_DIR="${LOG_DIR}"
    fi
    
    # Ensure LOG_FILE is absolute
    if [[ "$LOG_FILE" != /* ]]; then
        LOG_FILE="${LOG_DIR}/${LOG_FILE}"
    else
        LOG_FILE="${LOG_DIR}/$(basename "${LOG_FILE}")"
    fi
    
    # Set log level (can be overridden by environment variable)
    LOG_LEVEL="${LOG_LEVEL:-$DEFAULT_LOG_LEVEL}"
    
    # Store module name for use in log messages
    MODULE_NAME="$module_name"
    
    # Log initialization
    log_info "Logging system initialized for module: $module_name"
    log_info "Log type: $log_type"
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
