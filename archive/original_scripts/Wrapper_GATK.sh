#!/bin/bash

# =============================================================================
# INTERACTIVE WRAPPER FOR APPLE GATK PIPELINE - STAGE 1A
# =============================================================================
# AUTHORS: Phu Khang Ha
# DATE: $(date +%Y-%m-%d)
#
# DESCRIPTION:
# This interactive wrapper script automates the setup and execution of the
# Apple GATK Pipeline Stage 1A. It prompts the user for dataset information,
# automatically detects sample files, creates sample lists, and adjusts
# SLURM array parameters accordingly.
#
# FEATURES:
# - Interactive dataset name input with automatic path inference
# - Automatic sample counting and list generation
# - Optional user-provided sample list support
# - Dynamic SLURM array size adjustment
# - Comprehensive validation and error checking
# =============================================================================

# =============================================================================
# PROFESSIONAL LOGGING SYSTEM
# =============================================================================

# Colors for better user experience
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
PURPLE='\033[0;35m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

# Logging configuration
LOG_LEVELS=("DEBUG" "INFO" "WARN" "ERROR" "FATAL")
DEFAULT_LOG_LEVEL="INFO"
LOG_TO_CONSOLE=true
LOG_TO_FILE=true

# Initialize logging
init_logging() {
    local script_name=$(basename "$0" .sh)
    local timestamp=$(date '+%Y%m%d_%H%M%S')
    
    # Set log file path - separate from SLURM logs
    LOG_FILE="${script_name}_wrapper_${timestamp}.log"
    LOG_DIR="./logs/wrapper"
    
    # Create log directory if it doesn't exist
    mkdir -p "$LOG_DIR"
    LOG_FILE="$LOG_DIR/$LOG_FILE"
    
    # Set log level (can be overridden by environment variable)
    LOG_LEVEL="${LOG_LEVEL:-$DEFAULT_LOG_LEVEL}"
    
    # Log initialization with SLURM context awareness
    log_info "Interactive wrapper logging initialized"
    log_info "Wrapper log file: $LOG_FILE"
    log_info "Log level: $LOG_LEVEL"
    log_info "Note: SLURM job logs will be separate in /scratch/user/uqpha1/logs/"
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
    local level="$1"
    local message="$2"
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    local log_entry="[$timestamp] [$level] $message"
    
    # Check if we should log this level
    local current_level_num=$(get_log_level_num "$LOG_LEVEL")
    local message_level_num=$(get_log_level_num "$level")
    
    if [ $message_level_num -ge $current_level_num ]; then
        # Log to file
        if [ "$LOG_TO_FILE" = "true" ]; then
            echo "$log_entry" >> "$LOG_FILE"
        fi
        
        # Log to console with colors
        if [ "$LOG_TO_CONSOLE" = "true" ]; then
            case "$level" in
                "DEBUG") echo -e "${PURPLE}$log_entry${NC}";;
                "INFO")  echo -e "${BLUE}$log_entry${NC}";;
                "WARN")  echo -e "${YELLOW}$log_entry${NC}";;
                "ERROR") echo -e "${RED}$log_entry${NC}";;
                "FATAL") echo -e "${RED}$log_entry${NC}";;
                *)       echo "$log_entry";;
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

# SLURM-aware logging functions
log_slurm_submission() {
    local job_id="$1"
    local slurm_script="$2"
    local dataset_name="$3"
    local sample_count="$4"
    local slurm_log_dir="/scratch/user/uqpha1/logs/${dataset_name}"
    
    log_info "SLURM job submitted successfully"
    log_info "Job ID: $job_id"
    log_info "SLURM script: $slurm_script"
    log_info "Sample count: $sample_count"
    log_info "SLURM log directory: $slurm_log_dir"
    log_info "Monitor with: squeue -u \$USER"
    log_info "Cancel with: scancel $job_id"
    log_info "View SLURM logs: ls $slurm_log_dir/"
    
    # Create a job tracking file for easy reference
    local job_tracker="$LOG_DIR/job_${job_id}_${dataset_name}.info"
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

print_header() {
    echo ""
    echo "████████████████████████████████████████████████████████████████████████████████████████"
    echo "████████████████████████████████████████████████████████████████████████████████████████"
    echo "████████████████████████████████████████████████████████████████████████████████████████"
    echo ""
    echo "                            🚀 APPLE GATK PIPELINE - INTERACTIVE WRAPPER 🚀"
    echo ""
    echo "████████████████████████████████████████████████████████████████████████████████████████"
    echo "████████████████████████████████████████████████████████████████████████████████████████"
    echo "████████████████████████████████████████████████████████████████████████████████████████"
    echo ""
}

print_section() {
    local title=$1
    echo ""
    echo "┌─────────────────────────────────────────────────────────────────────────────────────┐"
    echo "│                           $title"
    echo "└─────────────────────────────────────────────────────────────────────────────────────┘"
    echo ""
}

# Function to validate directory structure
validate_directory_structure() {
    local base_path=$1
    local required_dirs=("1.FASTQ" "2.FASTQC_pre_trimmed" "3.FASTQC_post_trimmed" "4.BAM" "5.Individual_VCF" "6.Consolidated_GVCF")
    
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

# Function to check pipeline completion status for samples
check_pipeline_completion() {
    local rdm_base_path=$1
    local sample_list_file=$2
    
    log_info "Checking pipeline completion status..."
    print_message $BLUE "Checking pipeline completion status..."
    
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
            print_message $GREEN "✓ Completed: $sample"
        elif [ -f "$raw_gvcf" ] && [ -f "${raw_gvcf}.tbi" ]; then
            # Has GVCF but not genotyped VCF - partially complete
            incomplete_samples+=("$sample")
            log_warn "Sample partially complete: $sample (has GVCF, missing genotyped VCF)"
            print_message $YELLOW "⚠ Partially complete: $sample (has GVCF, missing genotyped VCF)"
        else
            incomplete_samples+=("$sample")
            log_debug "Sample incomplete: $sample"
            print_message $RED "✗ Incomplete: $sample"
        fi
    done < "$sample_list_file"
    
    print_message $CYAN "Pipeline Completion Summary:"
    print_message $GREEN "  Completed samples: ${#completed_samples[@]}"
    print_message $YELLOW "  Incomplete samples: ${#incomplete_samples[@]}"
    print_message $BLUE "  Total samples: $total_samples"
    
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

# Function to check Step 1B completion status
check_step1b_completion() {
    local rdm_base_path=$1
    
    log_info "Checking Step 1B (CombineGVCFs) completion status..."
    print_message $BLUE "Checking Step 1B completion status..."
    
    local consolidated_dir="${rdm_base_path}/6.Consolidated_GVCF"
    
    if [ ! -d "$consolidated_dir" ]; then
        log_warn "Step 1B directory does not exist: $consolidated_dir"
        print_message $RED "✗ Step 1B directory missing: $consolidated_dir"
        return 1
    fi
    
    # Check for consolidated VCF files (chromosome-based)
    local vcf_files=($(find "$consolidated_dir" -name "*.vcf.gz" -type f))
    local vcf_count=${#vcf_files[@]}
    
    if [ $vcf_count -eq 0 ]; then
        log_warn "No consolidated VCF files found in Step 1B directory"
        print_message $RED "✗ No consolidated VCF files found"
        return 1
    fi
    
    # Check if files have corresponding index files
    local complete_files=0
    for vcf_file in "${vcf_files[@]}"; do
        if [ -f "${vcf_file}.tbi" ]; then
            complete_files=$((complete_files + 1))
            local filename=$(basename "$vcf_file")
            print_message $GREEN "✓ Complete: $filename"
        else
            local filename=$(basename "$vcf_file")
            print_message $YELLOW "⚠ Missing index: $filename"
        fi
    done
    
    print_message $CYAN "Step 1B Completion Summary:"
    print_message $GREEN "  Complete VCF files: $complete_files"
    print_message $BLUE "  Total VCF files: $vcf_count"
    
    log_info "Step 1B completion summary - Complete: $complete_files, Total: $vcf_count"
    
    # Return completion status
    if [ $complete_files -eq $vcf_count ] && [ $vcf_count -gt 0 ]; then
        echo "complete"
    elif [ $complete_files -gt 0 ]; then
        echo "partial"
    else
        echo "incomplete"
    fi
}

# Function to run Step 1B (CombineGVCFs)
run_step1b() {
    local rdm_base_path=$1
    local dataset_name=$2
    
    log_info "Starting Step 1B (CombineGVCFs) execution..."
    print_message $BLUE "Starting Step 1B (CombineGVCFs)..."
    
    # Find the Step 1B script
    local step1b_script="Apple_GATK_Pipeline_1B_CombineGVCFs.sh"
    
    if [ ! -f "$step1b_script" ]; then
        log_error "Step 1B script not found: $step1b_script"
        print_message $RED "Error: Step 1B script not found: $step1b_script"
        print_message $YELLOW "Please ensure the script is in the current directory"
        return 1
    fi
    
    # Create SLURM submission script for Step 1B
    local slurm_script="Apple_GATK_1B_${dataset_name}_$(date +%Y%m%d_%H%M%S).sh"
    
    # SLURM configuration for Step 1B
    local job_name="Apple_GATK_1B_${dataset_name}"
    local account="a_qaafi_cas"
    local partition="general"
    local nodes=1
    local ntasks=1
    local cpus_per_task=6
    local time_limit="200:00:00"
    local memory="36G"
    local array_max=17  # Default for Step 1B (chromosome-based)
    
    # Ask user if they want to customize Step 1B SLURM parameters
    if confirm_action "Would you like to customize Step 1B SLURM parameters (default: 6 CPUs, 36G RAM, 200h time, array 0-17)?"; then
        print_message $CYAN "Step 1B SLURM Parameter Configuration:"
        
        cpus_per_task=$(get_user_input "Number of CPUs per task" "$cpus_per_task")
        memory=$(get_user_input "Memory allocation (e.g., 36G)" "$memory")
        time_limit=$(get_user_input "Time limit (format: HH:MM:SS)" "$time_limit")
        array_max=$(get_user_input "Array max (chromosome count)" "$array_max")
        
        print_message $GREEN "✓ Step 1B SLURM parameters configured:"
        print_message $YELLOW "  CPUs per task: $cpus_per_task"
        print_message $YELLOW "  Memory: $memory"
        print_message $YELLOW "  Time limit: $time_limit"
        print_message $YELLOW "  Array range: 0-${array_max}"
    fi
    
    # Create SLURM script for Step 1B
    cat > "$slurm_script" << EOF
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
#SBATCH -o /scratch/user/uqpha1/logs/${dataset_name}/1B_%A_%a_%x_%j.output
#SBATCH -e /scratch/user/uqpha1/logs/${dataset_name}/1B_%A_%a_%x_%j.error

# =============================================================================
# SLURM-GENERATED SCRIPT FOR APPLE GATK PIPELINE - STAGE 1B
# =============================================================================
# Generated by Interactive_GATK.sh on $(date)
# Dataset: ${dataset_name}
# Array Range: 0-${array_max}
# SLURM Config: ${cpus_per_task} CPUs, ${memory} RAM, ${time_limit} time
# =============================================================================

EOF
    
    # Append the Step 1B script content (skip the shebang line)
    tail -n +2 "$step1b_script" >> "$slurm_script"
    
    # Make it executable
    chmod +x "$slurm_script"
    
    print_message $GREEN "✓ Created Step 1B SLURM script: $slurm_script"
    print_message $CYAN "  Job name: $job_name"
    print_message $CYAN "  Array range: 0-${array_max}"
    print_message $CYAN "  SLURM config: ${cpus_per_task} CPUs, ${memory} RAM, ${time_limit} time"
    
    # Submit Step 1B job
    print_message $BLUE "Submitting Step 1B to SLURM..."
    
    # Create log directory for Step 1B
    mkdir -p "/scratch/user/uqpha1/logs/${dataset_name}"
    
    # Submit the job
    log_info "Submitting Step 1B SLURM job: $slurm_script with parameter: $rdm_base_path"
    local job_id=$(sbatch "$slurm_script" "$rdm_base_path" | awk '{print $4}')
    
    if [ $? -eq 0 ]; then
        # Use SLURM-aware logging for job submission
        log_slurm_submission "$job_id" "$slurm_script" "$dataset_name" "$((array_max + 1))"
        
        print_message $GREEN "✓ Step 1B submitted successfully!"
        print_message $CYAN "Job ID: $job_id"
        print_message $CYAN "Array jobs: 0-${array_max}"
        print_message $CYAN "SLURM Script: $slurm_script"
        print_message $YELLOW "Monitor with: squeue -u \$USER"
        print_message $YELLOW "Cancel with: scancel $job_id"
        print_message $YELLOW "View SLURM logs: ls /scratch/user/uqpha1/logs/${dataset_name}/1B_*"
        
        return 0
    else
        log_error "Failed to submit Step 1B to SLURM"
        print_message $RED "✗ Failed to submit Step 1B"
        return 1
    fi
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
    print_message $CYAN "Filtered sample list contents:"
    cat "$filtered_sample_list" | sed 's/^/  /'
    
    return ${#incomplete_samples[@]}
}

# Function to count samples and create sample list
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
    print_message $CYAN "Sample list contents:"
    cat "$sample_list_file" | sed 's/^/  /'
    
    return ${#valid_samples[@]}
}

# Function to create SLURM submission script
create_slurm_script() {
    local core_script=$1
    local sample_count=$2
    local dataset_name=$3
    local array_max=$((sample_count - 1))  # SLURM arrays are 0-indexed
    
    print_message $BLUE "Creating SLURM submission script..."
    
    # SLURM configuration options
    local job_name="Apple_GATK_1A_${dataset_name}"
    local account="a_qaafi_cas"
    local partition="general"
    local nodes=1
    local ntasks=1
    local cpus_per_task=10
    local time_limit="150:00:00"
    local memory="32G"
    
    # Ask user if they want to customize SLURM parameters
    if confirm_action "Would you like to customize SLURM parameters (default: 10 CPUs, 32G RAM, 150h time)?"; then
        print_message $CYAN "SLURM Parameter Configuration:"
        
        cpus_per_task=$(get_user_input "Number of CPUs per task" "$cpus_per_task")
        memory=$(get_user_input "Memory allocation (e.g., 32G)" "$memory")
        time_limit=$(get_user_input "Time limit (format: HH:MM:SS)" "$time_limit")
        
        print_message $GREEN "✓ SLURM parameters configured:"
        print_message $YELLOW "  CPUs per task: $cpus_per_task"
        print_message $YELLOW "  Memory: $memory"
        print_message $YELLOW "  Time limit: $time_limit"
    fi
    
    # Create a temporary script with SLURM directives
    local slurm_script="Apple_GATK_1A_${dataset_name}_$(date +%Y%m%d_%H%M%S).sh"
    
    cat > "$slurm_script" << EOF
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
#SBATCH -o /scratch/user/uqpha1/logs/${dataset_name}/%A_%a_%x_%j.output
#SBATCH -e /scratch/user/uqpha1/logs/${dataset_name}/%A_%a_%x_%j.error

# =============================================================================
# SLURM-GENERATED SCRIPT FOR APPLE GATK PIPELINE - STAGE 1A
# =============================================================================
# Generated by Interactive_1A.sh on $(date)
# Dataset: ${dataset_name}
# Sample Count: ${sample_count}
# Array Range: 0-${array_max}
# SLURM Config: ${cpus_per_task} CPUs, ${memory} RAM, ${time_limit} time
# =============================================================================

EOF
    
    # Append the core script content (skip the shebang line)
    tail -n +2 "$core_script" >> "$slurm_script"
    
    # Make it executable
    chmod +x "$slurm_script"
    
    print_message $GREEN "✓ Created SLURM script: $slurm_script"
    print_message $CYAN "  Job name: $job_name"
    print_message $CYAN "  Array range: 0-${array_max} (${sample_count} samples)"
    print_message $CYAN "  SLURM config: ${cpus_per_task} CPUs, ${memory} RAM, ${time_limit} time"
    print_message $CYAN "  Log directory: /scratch/user/uqpha1/logs/${dataset_name}/"
    
    echo "$slurm_script"
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

# Main execution
main() {
    # Initialize logging system first
    init_logging
    
    log_info "Starting Apple GATK Pipeline Interactive Wrapper"
    log_info "Script version: 2.0.0"
    log_info "User: $(whoami)"
    log_info "Host: $(hostname)"
    log_info "Working directory: $(pwd)"
    
    print_header
    
    print_section "📋 DATASET CONFIGURATION"
    
    # Get dataset name from user
    dataset_name=$(get_user_input "Enter dataset name (e.g., NCBI, WGS_Reference_panel)")
    
    log_info "User entered dataset name: $dataset_name"
    
    if [ -z "$dataset_name" ]; then
        log_error "Dataset name cannot be empty"
        print_message $RED "Error: Dataset name cannot be empty"
        exit 1
    fi
    
    # Infer RDM base path
    rdm_base_path="/QRISdata/Q8367/WGS_Reference_panel/${dataset_name}"
    
    log_info "Inferred RDM base path: $rdm_base_path"
    print_message $CYAN "Inferred RDM base path: $rdm_base_path"
    
    # Ask user to confirm or provide custom path
    custom_path=$(get_user_input "Confirm path or enter custom RDM base path" "$rdm_base_path")
    rdm_base_path="$custom_path"
    
    log_info "Final RDM base path: $rdm_base_path"
    
    # Validate the path exists
    if [ ! -d "$rdm_base_path" ]; then
        log_error "Path does not exist: $rdm_base_path"
        print_message $RED "Error: Path does not exist: $rdm_base_path"
        print_message $YELLOW "Please ensure the directory exists and contains the required structure"
        exit 1
    else
        log_info "Validating directory structure for: $rdm_base_path"
        validate_directory_structure "$rdm_base_path"
        
        # Check if critical directories are missing
        if [ ! -d "$rdm_base_path/1.FASTQ" ]; then
            log_error "Critical directory missing: $rdm_base_path/1.FASTQ"
            print_message $RED "Error: Critical directory missing: $rdm_base_path/1.FASTQ"
            print_message $YELLOW "This directory is required for the pipeline to run"
            exit 1
        fi
        log_info "Directory validation completed successfully"
    fi
    
    print_section "📁 SAMPLE LIST MANAGEMENT"
    
    # Ask user about sample list
    print_message $CYAN "Sample list options:"
    print_message $YELLOW "1. Auto-generate from FASTQ files in $rdm_base_path/1.FASTQ"
    print_message $YELLOW "2. Provide your own sample list file"
    
    while true; do
        read -p "Choose option (1/2): " option
        case $option in
            1)
                sample_list_mode="auto"
                break
                ;;
            2)
                sample_list_mode="manual"
                break
                ;;
            *)
                echo "Please choose 1 or 2"
                ;;
        esac
    done
    
    # Handle sample list creation
    sample_list_file="sample_list_${dataset_name}_$(date +%Y%m%d_%H%M%S).txt"
    
    if [ "$sample_list_mode" = "auto" ]; then
        print_message $BLUE "Auto-generating sample list..."
        sample_count=$(create_sample_list "$rdm_base_path/1.FASTQ" "$sample_list_file")
        if [ $? -ne 0 ]; then
            print_message $RED "Failed to create sample list. Exiting."
            exit 1
        fi
    else
        # Manual sample list
        while true; do
            sample_list_path=$(get_user_input "Enter path to your sample list file")
            if [ -f "$sample_list_path" ]; then
                cp "$sample_list_path" "$sample_list_file"
                sample_count=$(wc -l < "$sample_list_file")
                print_message $GREEN "✓ Using provided sample list with $sample_count samples"
                break
            else
                print_message $RED "File not found: $sample_list_path"
                if ! confirm_action "Try again?"; then
                    exit 1
                fi
            fi
        done
    fi
    
    print_section "🔍 PIPELINE COMPLETION CHECK"
    
    # Check if samples are already completed
    completion_status=$(check_pipeline_completion "$rdm_base_path" "$sample_list_file")
    
    case "$completion_status" in
        "all_complete")
            print_message $GREEN "🎉 All samples are already complete!"
            
            # Check Step 1B completion
            print_section "🔍 STEP 1B COMPLETION CHECK"
            step1b_status=$(check_step1b_completion "$rdm_base_path")
            
            case "$step1b_status" in
                "complete")
                    print_message $GREEN "🎉 Step 1B is also complete!"
                    if confirm_action "Would you like to rerun the entire pipeline anyway?"; then
                        print_message $YELLOW "Proceeding with full pipeline run..."
                    else
                        print_message $CYAN "Pipeline execution cancelled - all steps are complete"
                        exit 0
                    fi
                    ;;
                "partial"|"incomplete")
                    print_message $YELLOW "⚠ Step 1A is complete but Step 1B needs to be run"
                    if confirm_action "Would you like to run Step 1B (CombineGVCFs) now?"; then
                        print_message $BLUE "Running Step 1B..."
                        if run_step1b "$rdm_base_path" "$dataset_name"; then
                            print_message $GREEN "✓ Step 1B submitted successfully!"
                            print_message $CYAN "You can monitor the job with: squeue -u \$USER"
                            exit 0
                        else
                            print_message $RED "✗ Failed to submit Step 1B"
                            exit 1
                        fi
                    else
                        print_message $CYAN "Step 1B execution cancelled"
                        exit 0
                    fi
                    ;;
            esac
            ;;
        "partial_complete")
            print_message $YELLOW "⚠ Some samples are already complete"
            if confirm_action "Would you like to run only incomplete samples?"; then
                # Create filtered sample list
                filtered_sample_list="sample_list_${dataset_name}_incomplete_$(date +%Y%m%d_%H%M%S).txt"
                sample_count=$(create_filtered_sample_list "$rdm_base_path" "$sample_list_file" "$filtered_sample_list")
                if [ $? -ne 0 ]; then
                    print_message $GREEN "All samples are actually complete! Exiting."
                    exit 0
                fi
                sample_list_file="$filtered_sample_list"
                print_message $GREEN "Using filtered sample list with $sample_count incomplete samples"
            else
                print_message $YELLOW "Proceeding with full pipeline run (including completed samples)..."
            fi
            ;;
        "none_complete")
            print_message $BLUE "No samples are complete - proceeding with full pipeline run"
            ;;
    esac
    
    print_section "⚙️  SLURM CONFIGURATION"
    
    # Find the core script
    core_script="Apple_GATK_Pipeline_1a_CallVariantsPerSampleKHv2.sh"
    
    if [ ! -f "$core_script" ]; then
        print_message $RED "Error: Core script not found: $core_script"
        print_message $YELLOW "Please ensure the script is in the current directory"
        exit 1
    fi
    
    # Create SLURM submission script
    slurm_script=$(create_slurm_script "$core_script" "$sample_count" "$dataset_name")
    
    print_section "🚀 PIPELINE EXECUTION"
    
    # Display summary
    print_message $CYAN "Pipeline Configuration Summary:"
    print_message $YELLOW "  Dataset Name: $dataset_name"
    print_message $YELLOW "  RDM Base Path: $rdm_base_path"
    print_message $YELLOW "  Sample Count: $sample_count"
    print_message $YELLOW "  Sample List: $sample_list_file"
    print_message $YELLOW "  Core Script: $core_script"
    print_message $YELLOW "  SLURM Script: $slurm_script"
    print_message $YELLOW "  SLURM Array: 0-$((sample_count - 1))"
    
    echo ""
    
    # Confirm execution
    if confirm_action "Ready to submit the pipeline to SLURM?"; then
        log_info "User confirmed pipeline submission"
        print_message $GREEN "Submitting pipeline to SLURM..."
        
        # Create log directory
        mkdir -p "/scratch/user/uqpha1/logs/${dataset_name}"
        log_info "Created SLURM log directory: /scratch/user/uqpha1/logs/${dataset_name}"
        
        # Submit the job
        log_info "Submitting SLURM job: $slurm_script with parameters: $rdm_base_path $sample_list_file"
        job_id=$(sbatch "$slurm_script" "$rdm_base_path" "$sample_list_file" | awk '{print $4}')
        
        if [ $? -eq 0 ]; then
            # Use SLURM-aware logging for job submission
            log_slurm_submission "$job_id" "$slurm_script" "$dataset_name" "$sample_count"
            
            print_message $GREEN "✓ Pipeline submitted successfully!"
            print_message $CYAN "Job ID: $job_id"
            print_message $CYAN "Array jobs: 0-$((sample_count - 1))"
            print_message $CYAN "SLURM Script: $slurm_script"
            print_message $YELLOW "Monitor with: squeue -u \$USER"
            print_message $YELLOW "Cancel with: scancel $job_id"
            print_message $YELLOW "View SLURM logs: ls /scratch/user/uqpha1/logs/${dataset_name}/"
            print_message $YELLOW "View wrapper logs: $LOG_FILE"
            
            # Ask if user wants to automatically run Step 1B after Step 1A completes
            print_section "🔄 STEP 1B AUTOMATION"
            if confirm_action "Would you like to automatically run Step 1B (CombineGVCFs) after Step 1A completes?"; then
                print_message $BLUE "Step 1B will be automatically triggered after Step 1A completion"
                print_message $CYAN "You can monitor Step 1A progress with: squeue -u \$USER"
                print_message $CYAN "Step 1B will start automatically when Step 1A finishes"
                
                # Create a monitoring script for Step 1B automation
                local monitor_script="monitor_step1a_for_1b_${dataset_name}_${job_id}.sh"
                cat > "$monitor_script" << EOF
#!/bin/bash
# Monitor Step 1A completion and automatically run Step 1B
# Generated by Interactive_GATK.sh

echo "Monitoring Step 1A job: $job_id"
echo "Waiting for Step 1A to complete..."

# Wait for Step 1A to complete
while squeue -j $job_id >/dev/null 2>&1; do
    echo "Step 1A still running... (\$(date))"
    sleep 300  # Check every 5 minutes
done

echo "Step 1A completed! Checking Step 1B status..."

# Check if Step 1B is already complete
if [ -d "$rdm_base_path/6.Consolidated_GVCF" ]; then
    vcf_count=\$(find "$rdm_base_path/6.Consolidated_GVCF" -name "*.vcf.gz" -type f | wc -l)
    if [ \$vcf_count -gt 0 ]; then
        echo "Step 1B appears to be already complete (\$vcf_count VCF files found)"
        echo "Skipping Step 1B execution"
        exit 0
    fi
fi

echo "Step 1B not complete. Running Step 1B..."

# Run Step 1B
bash "$(pwd)/Interactive_GATK.sh" << 'INPUT'
$dataset_name
$rdm_base_path
INPUT

echo "Step 1B execution initiated"
EOF
                
                chmod +x "$monitor_script"
                print_message $GREEN "✓ Created monitoring script: $monitor_script"
                print_message $CYAN "To run monitoring in background: nohup bash $monitor_script &"
            else
                print_message $CYAN "Step 1B will not be automatically triggered"
                print_message $YELLOW "You can manually run Step 1B later using this script"
            fi
        else
            log_error "Failed to submit pipeline to SLURM"
            print_message $RED "✗ Failed to submit pipeline"
            exit 1
        fi
    else
        log_info "User cancelled pipeline submission"
        print_message $YELLOW "Pipeline submission cancelled"
        print_message $CYAN "You can run it manually with:"
        print_message $CYAN "sbatch $slurm_script $rdm_base_path $sample_list_file"
    fi
    
    log_info "Interactive wrapper completed successfully"
    print_message $GREEN "Interactive wrapper completed successfully!"
}

# Run main function
main "$@"
