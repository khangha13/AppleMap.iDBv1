#!/bin/bash -l
# =============================================================================
# SLURM DIRECTIVES REMOVED - NOW HANDLED BY INTERACTIVE_1A.SH
# =============================================================================
# This script is now designed to be submitted via sbatch by the Interactive_1A.sh wrapper
# The wrapper will add the appropriate SLURM directives based on the dataset and sample count

# =============================================================================
# APPLE GATK PIPELINE - STAGE 1A: PER-SAMPLE VARIANT CALLING
# =============================================================================
# AUTHORS: Phu Khang Ha
# DATE: 26/08/2025

#
# PIPELINE OVERVIEW:
# This script implements the GATK Best Practices workflow for variant calling
# in apple genomes. It processes each sample independently in parallel using
# SLURM array jobs.
# The calibration steps are based on the Howard et al knownsites 
#       
# =============================================================================
# DATA STORAGE STRUCTURE
# =============================================================================
# This pipeline uses an optimized storage strategy with minimal redundancy:
#
# 1. PERMANENT STORAGE (RDM - QRISdata):
#    /QRISdata/Q8367/WGS_Reference_panel/NCBI/
#    ├── 1.FASTQ                    # Raw sequencing data
#    ├── 2.FASTQC_pre_trimmed      # Quality control reports (before trimming)
#    ├── 3.FASTQC_post_trimmed     # Quality control reports (after trimming)
#    ├── 4.BAM                      # Aligned and processed BAM files
#    ├── 5.Individual_VCF          # Per-sample GVCF files
#    ├── 6.Consolidated_GVCF       # Combined GVCF files (Stage 1B)
#    └── 9.Metrics                  # Pipeline metrics files
#
# 2. TEMPORARY STORAGE (Scratch - Fast Processing - Essential Only):
#    /scratch/user/uqpha1/${today_folder}/
#    ├── 7.Backup                  # Step-by-step backups for resume functionality
#    └── 8.Time_logs               # Pipeline execution timing logs
#
# OPTIMIZATION: All outputs go directly to RDM permanent storage, eliminating
# redundant scratch copies and reducing I/O overhead significantly.
# =============================================================================


# PIPELINE OVERVIEW:
#
# Stage 1A (This script): Per-sample processing and GVCF generation
# Stage 1B (Separate script): Joint genotyping and variant filtering
# 
#
# WORKFLOW STEPS:
# 1. Raw read QC (FastQC)
# 2. Read trimming and quality filtering (Trimmomatic)
# 3. Alignment to reference + sort + index (BWA-MEM + Samtools)
# 4. Duplicate marking (GATK MarkDuplicates)
# 5. Base quality score recalibration (GATK BQSR)
# 6. Variant calling (GATK HaplotypeCaller with GVCF output)
# 7. Genotype GVCF to VCF conversion (GATK GenotypeGVCFs)
#
# INPUT: RDM base path containing standard folder structure + sample list
# OUTPUT: Per-sample GVCF and VCF files ready for analysis (saved directly to RDM)
# =============================================================================

# Get command line arguments
rdm_base_path=$1         # RDM base path containing the standard folder structure (e.g., /QRISdata/Q8367/NCBI)
sample_list=$2           # Text file containing list of sample IDs

# Extract dataset name from RDM path for scratch directory naming
today_dataset_folder=$(basename "$rdm_base_path")

# Define main working directory once (for all tasks)
readonly workdir="/scratch/user/uqpha1/${today_dataset_folder}_$(date +%Y%m)"

#Note
#The name in the sample list should match the name of the fastq files in the input data folder. This acts as a check point

# Validate command line arguments
if [ -z "$rdm_base_path" ] || [ -z "$sample_list" ]; then
    echo "Error: Missing required arguments"
    echo "Usage: sbatch script.sh <rdm_base_path> <sample_list>"
    echo "Example: sbatch script.sh /QRISdata/Q8367/NCBI sample_list.txt"
    echo ""
    echo "RDM base path should contain the standard folder structure:"
    echo "  /QRISdata/Q8367/NCBI/"
    echo "  ├── 1.FASTQ/"
    echo "  ├── 2.FASTQC_pre_trimmed/"
    echo "  ├── 3.FASTQC_post_trimmed/"
    echo "  ├── 4.BAM/"
    echo "  ├── 5.Individual_VCF/"
    echo "  └── 6.Consolidated_GVCF/"
    exit 1
fi

# Validate RDM base path exists and is a directory
if [ ! -d "$rdm_base_path" ]; then
    echo "Error: RDM base path '$rdm_base_path' does not exist or is not a directory"
    exit 1
fi

# Validate that RDM base path contains the expected folder structure
if [ ! -d "$rdm_base_path/1.FASTQ" ]; then
    echo "Error: RDM base path '$rdm_base_path' does not contain expected folder structure"
    echo "Expected: $rdm_base_path/1.FASTQ"
    echo "Please ensure the RDM path contains the standard folder structure"
    exit 1
fi

# Validate sample list file exists
# The script should be run from the home directory, therefore the sample list file should be in the home directory
if [ ! -f "$sample_list" ]; then
    echo "Error: Sample list file '$sample_list' does not exist"
    exit 1
fi

# Validate that sample list file is not empty
if [ ! -s "$sample_list" ]; then
    echo "Error: Sample list file '$sample_list' is empty"
    exit 1
fi

# =============================================================================
# PROFESSIONAL LOGGING SYSTEM FOR CORE PIPELINE
# =============================================================================
# This logging system works within SLURM's logging framework
# It provides structured logging for pipeline steps while respecting SLURM's output handling

# Logging configuration
PIPELINE_LOG_LEVELS=("DEBUG" "INFO" "WARN" "ERROR" "FATAL")
DEFAULT_PIPELINE_LOG_LEVEL="INFO"
PIPELINE_LOG_TO_CONSOLE=true
PIPELINE_LOG_TO_FILE=true

# Initialize pipeline logging
init_pipeline_logging() {
    local sample="$1"
    local timestamp=$(date '+%Y%m%d_%H%M%S')
    
    # Set pipeline log file path - use SLURM log directory structure
    # This ensures pipeline logs are co-located with SLURM job logs
    PIPELINE_LOG_FILE="pipeline_${sample}_${timestamp}.log"
    PIPELINE_LOG_DIR="/scratch/user/uqpha1/logs/${dataset_name}"
    
    # Create log directory if it doesn't exist (SLURM should have created it)
    mkdir -p "$PIPELINE_LOG_DIR"
    PIPELINE_LOG_FILE="$PIPELINE_LOG_DIR/$PIPELINE_LOG_FILE"
    
    # Set log level (can be overridden by environment variable)
    PIPELINE_LOG_LEVEL="${PIPELINE_LOG_LEVEL:-$DEFAULT_PIPELINE_LOG_LEVEL}"
    
    # Log initialization
    pipeline_log_info "Pipeline logging initialized for sample: $sample"
    pipeline_log_info "Pipeline log file: $PIPELINE_LOG_FILE"
    pipeline_log_info "Log level: $PIPELINE_LOG_LEVEL"
    pipeline_log_info "SLURM Array Task ID: $SLURM_ARRAY_TASK_ID"
    pipeline_log_info "SLURM Job ID: $SLURM_JOB_ID"
    pipeline_log_info "Dataset: $dataset_name"
}

# Get numeric log level
get_pipeline_log_level_num() {
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

# Core pipeline logging function
pipeline_log_message() {
    local level="$1"
    local message="$2"
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    local log_entry="[$timestamp] [$level] [Sample:$SAMPLE] [Task:$SLURM_ARRAY_TASK_ID] $message"
    
    # Check if we should log this level
    local current_level_num=$(get_pipeline_log_level_num "$PIPELINE_LOG_LEVEL")
    local message_level_num=$(get_pipeline_log_level_num "$level")
    
    if [ $message_level_num -ge $current_level_num ]; then
        # Log to file
        if [ "$PIPELINE_LOG_TO_FILE" = "true" ]; then
            echo "$log_entry" >> "$PIPELINE_LOG_FILE"
        fi
        
        # Log to console (SLURM will capture this)
        if [ "$PIPELINE_LOG_TO_CONSOLE" = "true" ]; then
            echo "$log_entry"
        fi
    fi
}

# Convenience pipeline logging functions
pipeline_log_debug() { pipeline_log_message "DEBUG" "$1"; }
pipeline_log_info()  { pipeline_log_message "INFO"  "$1"; }
pipeline_log_warn()  { pipeline_log_message "WARN"  "$1"; }
pipeline_log_error() { pipeline_log_message "ERROR" "$1"; }
pipeline_log_fatal() { pipeline_log_message "FATAL" "$1"; }

# Legacy logging function (for backward compatibility)
log_message() {
    local message="$1"
    pipeline_log_info "$message"
    
    # Also provide the old formatted output for visual consistency
    echo ""
    echo "=================================================================================="
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $message"
    echo ""
}

# Function to setup shared reference files (only run by task 0)
setup_shared_reference_files() {
    if [ "$SLURM_ARRAY_TASK_ID" -eq 0 ]; then
        log_message "Setting up shared reference files for all array jobs"
        
        # Check if reference files already exist
        if [ ! -f "$SHARED_REF_GENOME_DIR/$REF_GENOME" ]; then
            log_message "Copying reference genome and indices to shared directory"
            if ! rsync -rhivPt "${REF_GENOME_DIR}/" "$SHARED_REF_GENOME_DIR/"; then
                error_exit "Failed to copy reference genome files to shared directory"
            fi
            
            log_message "Copying known sites to shared directory"
            if ! rsync -rhivPt "${KNOWN_SITES_DIR}/" "$SHARED_KNOWN_SITES_DIR/"; then
                error_exit "Failed to copy known sites files to shared directory"
            fi
            
            log_message "Copying adapter file to shared directory"
            if ! rsync -rhivPt "${ADAPTER_FILE}" "$SHARED_ADAPTER_DIR/"; then
                error_exit "Failed to copy adapter file to shared directory"
            fi
            
            log_message "Shared reference files setup completed"
        else
            log_message "Shared reference files already exist, skipping copy"
        fi
    fi
}

# Error handling function - exits with error message (defined once)

# =============================================================================
# INITIALISATION PHASE - DIRECTORY CREATION
# =============================================================================
# Create all required directories for the pipeline
# This ensures the pipeline has proper storage locations before processing

# Only create directories if this is the first array task (SLURM_ARRAY_TASK_ID = 0)
if [ "$SLURM_ARRAY_TASK_ID" -eq 0 ]; then
    log_message "Initializing pipeline directories for run: $today_dataset_folder"
    
    # Function to check and create directory if it doesn't exist
    create_dir_if_not_exists() {
        local dir="$1"
        local dir_type="$2"
        
        if [ -d "$dir" ]; then
            log_message "Directory already exists: $dir_type ($dir)"
        else
            log_message "Creating directory: $dir_type ($dir)"
            mkdir -p "$dir"
            if [ $? -eq 0 ]; then
                log_message "Successfully created: $dir_type"
            else
                error_exit "Failed to create directory: $dir_type ($dir)"
            fi
        fi
    }
    
    # Create main directories with existence checks
    create_dir_if_not_exists "$workdir" "Scratch working directory"
    
    # Create subdirectories in scratch with existence checks (only essential directories)
    create_dir_if_not_exists "$workdir/7.Backup" "Step-by-step backups directory"
    create_dir_if_not_exists "$workdir/8.Time_logs" "Timing logs directory"
    
    # Create shared reference directories for all array jobs
    create_dir_if_not_exists "$workdir/9.Reference_genome_dir" "Shared reference genome directory"
    create_dir_if_not_exists "$workdir/10.Knownsites_dir" "Shared known sites directory"
    create_dir_if_not_exists "$workdir/11.Adapter_file" "Shared adapter file directory"
    
    log_message "Created directory structure:"
    log_message "  Scratch: $workdir"
    log_message "  RDM base: $rdm_base_path"
    log_message "Initialization completed successfully"
    
    # Setup shared reference files (only done by task 0)
    setup_shared_reference_files
fi

# Wait for initialization to complete (only relevant for array jobs)
if [ "$SLURM_ARRAY_TASK_ID" -gt 0 ]; then
    # Small delay to ensure directories and reference files are created by task 0
    sleep 10
fi

# =============================================================================
# MODULE LOADING SECTION
# =============================================================================
# Load all required bioinformatics tools and their specific versions
# These versions have been tested and are compatible with each other

module load fastqc/0.11.9-java-11              # Quality control for raw reads
module load trimmomatic/0.39-java-11           # Read trimming and adapter removal
module load bwa/0.7.17-gcccore-11.3.0          # Burrows-Wheeler Aligner for read mapping
module load samtools/1.16.1-gcc-11.3.0         # SAM/BAM file manipulation and indexing
module load gatk/4.3.0.0-gcccore-11.3.0-java-11 # Genome Analysis Toolkit for variant calling
module load picard                              # Java tools for working with sequencing data



# =============================================================================
# PIPELINE PARAMETERS AND CONSTANTS
# =============================================================================
# These parameters have been optimized for apple genome analysis
# Adjust these values based on your specific dataset and requirements

readonly NUM_THREADS=12                         # Number of CPU threads for parallel processing
readonly MEMORY="24G"                          # Memory allocation for Java-based tools
readonly MIN_QUALITY=20                        # Minimum base quality score for trimming
readonly SLIDING_WINDOW="3:15"                 # Sliding window size:quality threshold for trimming
readonly MIN_LENGTH=36                         # Minimum read length after trimming
readonly TRIM_LEADING=20                       # Trim low-quality bases from read start
readonly TRIM_TRAILING=20                      # Trim low-quality bases from read end
readonly ILLUMINA_CLIP_SETTINGS="2:30:3:1:True" # Adapter clipping parameters (mismatches:palindrome:simple:minAdapterLength:keepBothReads)

# =============================================================================
# REFERENCE FILES AND PATHS
# =============================================================================
# Define all reference files and their locations
# These files are specific to the apple genome analysis

readonly REF_GENOME_DIR="/QRISdata/Q8367/Reference_Genome"   # Directory containing reference genome and indices
readonly REF_GENOME="GDDH13_1-1_formatted.fasta"                    # Apple reference genome (GDDH13 v1.1)
readonly KNOWN_SITES_DIR="/QRISdata/Q8367/Known_Sites"       # Directory containing known variant sites
readonly KNOWN_SITES="Final_DB_Known-sites_Filtered_by_480K_Unique_HD.vcf"  # Known variant sites for BQSR
readonly ADAPTER_FILE="/QRISdata/Q8367/Liao_2021_adapter/TruSeq2-PE.fa"  # Illumina adapter sequences for Liao et al. 2021




# =============================================================================
# SAMPLE PROCESSING SETUP
# =============================================================================
# Each SLURM array task processes one sample from the sample list
# Array IDs start from 0, so we add 1 to get the correct line number

SAMPLE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$sample_list")

if [ -z "$SAMPLE" ]; then
    echo "Error: No sample found for task ID ${SLURM_ARRAY_TASK_ID} (line $((SLURM_ARRAY_TASK_ID + 1)))"
    exit 1
fi

# Extract dataset name from RDM base path for logging
dataset_name=$(basename "$RDM_BASE")

# Initialize pipeline logging for this sample
init_pipeline_logging "$SAMPLE"

# =============================================================================
# WORKING DIRECTORY STRUCTURE
# =============================================================================
# Use the directories created during initialization
# Each sample gets its own subdirectories for different processing stages

# Define main directories (already created in initialization)

# Shared reference directories for all array jobs
readonly SHARED_REF_GENOME_DIR="${workdir}/9.Reference_genome_dir"
readonly SHARED_KNOWN_SITES_DIR="${workdir}/10.Knownsites_dir"
readonly SHARED_ADAPTER_DIR="${workdir}/11.Adapter_file"

# Define all directory paths for different pipeline stages (RDM-only storage)
readonly BACKUP_DIR="${workdir}/7.Backup/${SAMPLE}"                   # Step-by-step backups. Each sample will have its own backup directory.

# RDM (Permanent) storage paths - using the input RDM base path
readonly RDM_BASE="$rdm_base_path"
readonly RDM_FASTQ="${RDM_BASE}/1.FASTQ"
readonly RDM_QC_PRE="${RDM_BASE}/2.FASTQC_pre_trimmed"
readonly RDM_QC_POST="${RDM_BASE}/3.FASTQC_post_trimmed"
readonly RDM_BAM="${RDM_BASE}/4.BAM"
readonly RDM_VCF="${RDM_BASE}/5.Individual_VCF"
readonly RDM_CONSOLIDATED="${RDM_BASE}/6.Consolidated_GVCF"
readonly RDM_METRICS="${RDM_BASE}/9.Metrics"

# =============================================================================
# READ GROUP INFORMATION
# =============================================================================
# Read groups are essential for GATK variant calling
# They help identify the source of each read and enable proper variant calling

readonly READ_GROUP="@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:LIB1\tPU:UNIT1"
# ID: Unique identifier for this read group
# SM: Sample name
# PL: Platform (ILLUMINA)
# LB: Library identifier
# PU: Platform unit (sequencing run identifier)

# =============================================================================
# SAMPLE-SPECIFIC DIRECTORY CREATION
# =============================================================================
# Create sample-specific subdirectories within the initialized structure
# These are created per sample during processing

# Display path information for user reference
log_message "Processing sample: $SAMPLE"
log_message "RDM base path: $rdm_base_path"
log_message "Sample data location: ${RDM_FASTQ}/${SAMPLE}"
log_message "Working directory: $workdir"

echo ""
echo ""
echo "████████████████████████████████████████████████████████████████████████████████████████"
echo "████████████████████████████████████████████████████████████████████████████████████████"
echo "████████████████████████████████████████████████████████████████████████████████████████"
echo ""
echo "                                    🚀 APPLE GATK PIPELINE STARTING 🚀"
echo ""
echo "████████████████████████████████████████████████████████████████████████████████████████"
echo "████████████████████████████████████████████████████████████████████████████████████████"
echo "████████████████████████████████████████████████████████████████████████████████████████"
echo ""
    echo "    🧬 Sample: $SAMPLE"
    echo "    📁 RDM Base Path: $rdm_base_path"
    echo "    🏠 Working Directory: $workdir"
    echo "    🕐 Started: $(date)"
echo ""
echo "████████████████████████████████████████████████████████████████████████████████████████"
echo "████████████████████████████████████████████████████████████████████████████████████████"
echo "████████████████████████████████████████████████████████████████████████████████████████"
echo ""
echo ""

# Change to working directory (TMPDIR for fast local processing)
cd $TMPDIR || error_exit "Failed to change to TMPDIR"

# Create sample-specific directories (RDM-only storage)
mkdir -p "$BACKUP_DIR"
mkdir -p "$RDM_BAM/${SAMPLE}" "$RDM_QC_PRE/${SAMPLE}" "$RDM_QC_POST/${SAMPLE}" \
         "$RDM_VCF" "$RDM_METRICS/${SAMPLE}"


# =============================================================================
# TIMING FUNCTIONS AND ERROR HANDLING
# =============================================================================
# Functions to track execution time for each step and overall pipeline
# Includes error handling to report timing even when script fails

# Global timing variables
PIPELINE_START_TIME=$(date +%s)
STEP_START_TIME=0
STEP_END_TIME=0
CURRENT_STEP_NAME=""
STEP_TIMES_FILE="${workdir}/8.Time_logs/step_times_${SAMPLE}.log"

# Function to log step timing to file for persistence
log_step_time() {
    local step_name="$1"
    local timestamp="$2"
    local status="$3"
    echo "[$timestamp] $status: $step_name" >> "$STEP_TIMES_FILE"
}

# Function to report total time on script exit (success or failure)
report_total_time_on_exit() {
    local exit_code=$?
    local pipeline_end_time=$(date +%s)
    local total_duration=$((pipeline_end_time - PIPELINE_START_TIME))
    local hours=$((total_duration / 3600))
    local minutes=$(((total_duration % 3600) / 60))
    local seconds=$((total_duration % 60))
    
    echo ""
    echo ""
    echo "████████████████████████████████████████████████████████████████████████████████████████"
    echo "████████████████████████████████████████████████████████████████████████████████████████"
    echo "████████████████████████████████████████████████████████████████████████████████████████"
    echo ""
    if [ $exit_code -eq 0 ]; then
        echo "                                    🎉 PIPELINE COMPLETION SUMMARY 🎉"
    else
        echo "                                    ⚠️  PIPELINE ERROR SUMMARY ⚠️"
    fi
    echo ""
    echo "████████████████████████████████████████████████████████████████████████████████████████"
    echo "████████████████████████████████████████████████████████████████████████████████████████"
    echo "████████████████████████████████████████████████████████████████████████████████████████"
    echo ""
    echo "    📊 Sample: $SAMPLE"
    echo "    🕐 Pipeline Started:  $(date -d @$PIPELINE_START_TIME)"
    echo "    🕐 Pipeline Ended:    $(date -d @$pipeline_end_time)"
    echo ""
    
    if [ $hours -gt 0 ]; then
        echo "    ⏱️  Total Pipeline Time: ${hours}h ${minutes}m ${seconds}s"
    elif [ $minutes -gt 0 ]; then
        echo "    ⏱️  Total Pipeline Time: ${minutes}m ${seconds}s"
    else
        echo "    ⏱️  Total Pipeline Time: ${seconds}s"
    fi
    
    if [ $exit_code -ne 0 ]; then
        echo "    ❌ Exit Code: $exit_code"
        echo "    🔍 Check logs for error details"
    fi
    
    echo ""
    echo "    📁 Step Times Log: $STEP_TIMES_FILE"
    echo "    📁 Working Directory: $workdir"
    echo ""
    echo "████████████████████████████████████████████████████████████████████████████████████████"
    echo "████████████████████████████████████████████████████████████████████████████████████████"
    echo "████████████████████████████████████████████████████████████████████████████████████████"
    echo ""
    echo ""
    
    # Log the final timing to file
    log_step_time "PIPELINE_COMPLETE" "$(date '+%Y-%m-%d %H:%M:%S')" "EXIT_CODE_${exit_code}"
}

# Set up trap to call timing function on script exit (success or failure)
trap report_total_time_on_exit EXIT

# Function to start timing a step with enhanced formatting
start_step_timer() {
    local step_name="$1"
    STEP_START_TIME=$(date +%s)
    CURRENT_STEP_NAME="$step_name"
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    
    # Log step start with professional logging
    pipeline_log_info "Starting step: $step_name"
    
    echo ""
    echo "╔══════════════════════════════════════════════════════════════════════════════════╗"
    echo "║                           STARTING: $step_name"
    echo "║                           Started at: $timestamp"
    echo "╚══════════════════════════════════════════════════════════════════════════════════╝"
    echo ""
    
    # Log step start to file
    log_step_time "$step_name" "$timestamp" "STARTED"
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
        pipeline_log_info "Completed step: $step_name (Duration: ${hours}h ${minutes}m ${seconds}s)"
    elif [ $minutes -gt 0 ]; then
        pipeline_log_info "Completed step: $step_name (Duration: ${minutes}m ${seconds}s)"
    else
        pipeline_log_info "Completed step: $step_name (Duration: ${seconds}s)"
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
    
    # Log step completion to file
    log_step_time "$step_name" "$timestamp" "COMPLETED_${hours}h_${minutes}m_${seconds}s"
}

# Function to calculate and display total pipeline time with enhanced formatting
show_total_time() {
    local pipeline_end_time=$(date +%s)
    local total_duration=$((pipeline_end_time - PIPELINE_START_TIME))
    local hours=$((total_duration / 3600))
    local minutes=$(((total_duration % 3600) / 60))
    local seconds=$((total_duration % 60))
    
    echo ""
    echo ""
    echo "████████████████████████████████████████████████████████████████████████████████████████"
    echo "████████████████████████████████████████████████████████████████████████████████████████"
    echo "████████████████████████████████████████████████████████████████████████████████████████"
    echo ""
    echo "                                    🎉 PIPELINE COMPLETION SUMMARY 🎉"
    echo ""
    echo "████████████████████████████████████████████████████████████████████████████████████████"
    echo "████████████████████████████████████████████████████████████████████████████████████████"
    echo "████████████████████████████████████████████████████████████████████████████████████████"
    echo ""
    echo "    📊 Sample Processed: $SAMPLE"
    echo ""
    echo "    🕐 Pipeline Started:  $(date -d @$PIPELINE_START_TIME)"
    echo "    🕐 Pipeline Completed: $(date -d @$pipeline_end_time)"
    echo ""
    
    if [ $hours -gt 0 ]; then
        echo "    ⏱️  Total Pipeline Time: ${hours}h ${minutes}m ${seconds}s"
    elif [ $minutes -gt 0 ]; then
        echo "    ⏱️  Total Pipeline Time: ${minutes}m ${seconds}s"
    else
        echo "    ⏱️  Total Pipeline Time: ${seconds}s"
    fi
    echo ""
    echo "    📁 GVCF File Location: $RDM_VCF/${SAMPLE}_raw.g.vcf.gz"
    echo "    📁 Genotyped VCF Location: $RDM_VCF/${SAMPLE}_genotyped.vcf.gz"
    echo "    📁 Analysis BAM Location: $RDM_BAM/${SAMPLE}/${SAMPLE}_recal.bam"
    echo "    📁 Backup Location: $BACKUP_DIR"
    echo "    📁 RDM Permanent Storage: $RDM_BASE"
    echo ""
    echo "████████████████████████████████████████████████████████████████████████████████████████"
    echo "████████████████████████████████████████████████████████████████████████████████████████"
    echo "████████████████████████████████████████████████████████████████████████████████████████"
    echo ""
    echo ""
}

# =============================================================================

# Error handling function - exits with error message and timing info
error_exit() {
    local error_message="$1"
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    local current_time=$(date +%s)
    local runtime=$((current_time - PIPELINE_START_TIME))
    local hours=$((runtime / 3600))
    local minutes=$(((runtime % 3600) / 60))
    local seconds=$((runtime % 60))
    
    # Log error with professional logging
    pipeline_log_fatal "Pipeline error: $error_message"
    if [ -n "$CURRENT_STEP_NAME" ]; then
        pipeline_log_fatal "Failed at step: $CURRENT_STEP_NAME"
    fi
    pipeline_log_error "Runtime before error: ${hours}h ${minutes}m ${seconds}s"
    
    echo ""
    echo "████████████████████████████████████████████████████████████████████████████████████████"
    echo "████████████████████████████████████████████████████████████████████████████████████████"
    echo "████████████████████████████████████████████████████████████████████████████████████████"
    echo ""
    echo "                                    ❌ PIPELINE ERROR OCCURRED ❌"
    echo ""
    echo "████████████████████████████████████████████████████████████████████████████████████████"
    echo "████████████████████████████████████████████████████████████████████████████████████████"
    echo "████████████████████████████████████████████████████████████████████████████████████████"
    echo ""
    echo "    📊 Sample: $SAMPLE"
    echo "    🕐 Error Time: $timestamp"
    echo "    🔍 Error Message: $error_message"
    echo ""
    if [ -n "$CURRENT_STEP_NAME" ]; then
        echo "    ⚠️  Failed at Step: $CURRENT_STEP_NAME"
    fi
    echo ""
    if [ $hours -gt 0 ]; then
        echo "    ⏱️  Runtime Before Error: ${hours}h ${minutes}m ${seconds}s"
    elif [ $minutes -gt 0 ]; then
        echo "    ⏱️  Runtime Before Error: ${minutes}m ${seconds}s"
    else
        echo "    ⏱️  Runtime Before Error: ${seconds}s"
    fi
    echo ""
    echo "    📁 Step Times Log: $STEP_TIMES_FILE"
    echo "    📁 Pipeline Log: $PIPELINE_LOG_FILE"
    echo "    📁 Working Directory: $workdir"
    echo ""
    echo "████████████████████████████████████████████████████████████████████████████████████████"
    echo "████████████████████████████████████████████████████████████████████████████████████████"
    echo "████████████████████████████████████████████████████████████████████████████████████████"
    echo ""
    echo ""
    
    # Log error to timing file
    log_step_time "ERROR_OCCURRED" "$timestamp" "ERROR_${hours}h_${minutes}m_${seconds}s"
    
    echo "Error: $error_message" >&2
    exit 1
}

# =============================================================================
# SIMPLIFIED BACKUP SYSTEM (OPTIMIZED)
# =============================================================================
# Each sample has a single backup folder that mirrors TMPDIR after each step
# RDM-based detection determines the current pipeline step
# ECONOMIC OPTIMIZATION: Backup only starts from Step 3 (BWA alignment)
# Steps 1-2 (QC, trimming) are not backed up as they can be regenerated from RDM
# OPTIMIZED BACKUP: Only sample-specific files are backed up (reference files excluded)
# Reference files are now stored in shared directories to avoid duplication

# Function to detect current step based on RDM storage files
detect_current_step() {
    local sample="$SAMPLE"
    
    # Helper function to check if files exist in RDM
    check_rdm_files_exist() {
        for pattern in "$@"; do
            if ! ls $pattern >/dev/null 2>&1; then
                return 1
            fi
        done
        return 0
    }
    
    # Check steps in reverse order (most complete first)
    if check_rdm_files_exist "${RDM_VCF}/${sample}_genotyped.vcf.gz"; then
        echo "7"  # GenotypeGVCFs completed (Step 7)
    elif check_rdm_files_exist "${RDM_VCF}/${sample}_raw.g.vcf.gz"; then
        echo "6"  # HaplotypeCaller completed (Step 6)
    elif check_rdm_files_exist "${RDM_BAM}/${sample}/${sample}_recal.bam"; then
        echo "5"  # BQSR completed
    elif check_rdm_files_exist "${RDM_BAM}/${sample}/${sample}_dedup.bam"; then
        echo "4"  # MarkDuplicates completed
    elif check_rdm_files_exist "${RDM_BAM}/${sample}/${sample}_sorted.bam"; then
        echo "3"  # BWA alignment completed
    elif check_rdm_files_exist "${RDM_QC_POST}/${sample}"/*.html "${RDM_QC_POST}/${sample}"/*.zip; then
        echo "2"  # Post-trimming QC completed (Step 2 finished)
    elif check_rdm_files_exist "${RDM_QC_PRE}/${sample}"/*.html "${RDM_QC_PRE}/${sample}"/*.zip; then
        echo "1"  # Pre-trimming QC completed (Step 1 finished)
    else
        echo "0"  # No results found, start from step 1
    fi
}

# Simplified backup function - mirrors TMPDIR to backup folder (comprehensive backup)
backup_current_state() {
    local step_name="$1"
    local backup_dir="$BACKUP_DIR"
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    
    log_message "Creating comprehensive backup for step: $step_name"
    
    # Backup only sample-specific files, exclude reference files (now in shared directories)
    # This includes: BAM files, VCF files, FASTQ files, etc.
    if ls * 2>/dev/null; then
        # Exclude reference files from backup since they're now in shared directories
        rsync -rhivPt --exclude="*.fasta*" --exclude="*.vcf*" --exclude="*.fa" * "$backup_dir/" 2>/dev/null || true
        log_message "Backup completed for step: $step_name (reference files excluded)"
        
        # Log backup completion to timing file
        log_step_time "BACKUP_${step_name}" "$timestamp" "BACKUP_COMPLETED"
    else
        log_message "No files found to backup for step: $step_name"
        log_step_time "BACKUP_${step_name}" "$timestamp" "NO_FILES_TO_BACKUP"
    fi
}

# Function to restore files from backup
restore_from_backup() {
    local backup_dir="$BACKUP_DIR"
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    
    if [ -d "$backup_dir" ] && [ "$(ls -A "$backup_dir" 2>/dev/null)" ]; then
        log_message "Restoring files from backup directory: $backup_dir"
        log_message "Files to restore:"
        ls -la "$backup_dir" || log_message "Could not list backup files"
        
        rsync -rhivPt "$backup_dir"/* . 2>/dev/null || true
        log_message "Restoration completed"
        
        # Verify restoration by listing current directory
        log_message "Files in TMPDIR after restoration:"
        ls -l
        
        # Log restoration to timing file
        log_step_time "RESTORE_FROM_BACKUP" "$timestamp" "RESTORATION_COMPLETED"
    else
        log_message "No backup found to restore from"
        log_step_time "RESTORE_FROM_BACKUP" "$timestamp" "NO_BACKUP_FOUND"
    fi
}

# Function to prune obsolete files from backup directory to save space
prune_backup_files() {
    local description="$1"
    shift
    local files_to_remove=("$@")
    if [ ${#files_to_remove[@]} -gt 0 ]; then
        log_message "Pruning obsolete files from backup for: $description"
        for f in "${files_to_remove[@]}"; do
            rm -f "$BACKUP_DIR/$f" 2>/dev/null || true
        done
    fi
}

# =============================================================================
# STEP CHECKING AND RESUME FUNCTIONALITY (RDM-BASED)
# =============================================================================
# Use RDM-based detection to determine current step and resume from there

# Initialize timing log file
log_message "Initializing timing log: $STEP_TIMES_FILE"
log_step_time "PIPELINE_START" "$(date '+%Y-%m-%d %H:%M:%S')" "PIPELINE_INITIALIZED"

# Check completion status using RDM file detection
log_message "Checking pipeline completion status using RDM file detection..."

# Detect current step from RDM files
CURRENT_STEP=$(detect_current_step)
log_message "Detected current step: $CURRENT_STEP"

# Debug: Show what files were found in RDM directories
log_message "Checking RDM storage for existing files..."

# Check each RDM directory for sample files
for dir in "${RDM_QC_PRE}/${SAMPLE}" "${RDM_QC_POST}/${SAMPLE}" "${RDM_BAM}/${SAMPLE}" "${RDM_VCF}"; do
    if [ -d "$dir" ]; then
        log_message "Files in $dir:"
        ls -la "$dir" | grep "$SAMPLE" || log_message "No sample files found"
    else
        log_message "Directory does not exist: $dir"
    fi
done

# Determine starting step based on RDM file detection
if [ "$CURRENT_STEP" -eq 7 ]; then
    log_message "All steps completed! Pipeline finished for sample: $SAMPLE"
    exit 0
elif [ "$CURRENT_STEP" -eq 6 ]; then
    STARTING_STEP=7
    log_message "Resuming from Step 7: GenotypeGVCFs"
elif [ "$CURRENT_STEP" -eq 5 ]; then
    STARTING_STEP=6
    log_message "Resuming from Step 6: Variant Calling"
elif [ "$CURRENT_STEP" -eq 4 ]; then
    STARTING_STEP=5
    log_message "Resuming from Step 5: Base Quality Score Recalibration"
elif [ "$CURRENT_STEP" -eq 3 ]; then
    STARTING_STEP=4
    log_message "Resuming from Step 4: Duplicate Marking"
elif [ "$CURRENT_STEP" -eq 2 ]; then
    STARTING_STEP=3
    log_message "Resuming from Step 3: BWA Alignment (QC steps already completed)"
elif [ "$CURRENT_STEP" -eq 1 ]; then
    STARTING_STEP=2
    log_message "Resuming from Step 2: Read Trimming (pre-trimming QC already completed)"
else
    STARTING_STEP=1
    log_message "Starting from Step 1: Quality Control (no previous results found)"
fi

# Debug: Confirm the starting step
log_message "Starting step determined: $STARTING_STEP"

# =============================================================================
# FILE PREPARATION PHASE
# =============================================================================
# Copy all necessary reference files and data to TMPDIR for fast processing
# TMPDIR is typically on fast local storage for better I/O performance


# Copy reference files from shared directories to TMPDIR (Steps 1-2, 7)
if [ "$STARTING_STEP" -le 2 ] || [ "$STARTING_STEP" -eq 7 ]; then
    log_message "Copying reference files from shared directories to \$TMPDIR (needed for Steps 1-2 and Step 7)"
    
    # Create symlinks to shared reference files for fast access
    log_message "Creating symlinks to shared reference files"
    ln -sf "$SHARED_REF_GENOME_DIR"/* .
    ln -sf "$SHARED_KNOWN_SITES_DIR"/* .
    ln -sf "$SHARED_ADAPTER_DIR"/* .
    
    # Update paths for GATK commands
    readonly KNOWN_SITES_PATH="$SHARED_KNOWN_SITES_DIR/${KNOWN_SITES}"
    readonly ADAPTER_PATH="$SHARED_ADAPTER_DIR/$(basename "$ADAPTER_FILE")"
    
    log_message "Reference files linked successfully"
else
    log_message "Skipping reference files copy - will be restored from backup (Steps 3-6)"
fi

# Copy raw reads from RDM to $TMPDIR only if needed (Steps 1-2)
if [ "$STARTING_STEP" -le 2 ]; then
    log_message "Copying raw reads from RDM to \$TMPDIR (needed for Steps 1-2)"
    
    # Check if paired-end sample files exist
    if [ ! -f "${RDM_FASTQ}/${SAMPLE}_1.fastq.gz" ] || [ ! -f "${RDM_FASTQ}/${SAMPLE}_2.fastq.gz" ]; then
        error_exit "Paired-end sample files ${RDM_FASTQ}/${SAMPLE}_1.fastq.gz and ${RDM_FASTQ}/${SAMPLE}_2.fastq.gz do not exist in RDM"
    fi
    
    if ! rsync -rhivPt "${RDM_FASTQ}/${SAMPLE}_1.fastq.gz" "${RDM_FASTQ}/${SAMPLE}_2.fastq.gz" .; then
        error_exit "Failed to copy paired-end raw reads from RDM using rsync"
    fi
    
    # Keep the paired-end fastq.gz files compressed for downstream processing
    log_message "Keeping paired-end raw reads compressed (.gz format)"
else
    log_message "Skipping raw reads copy - not needed for starting step $STARTING_STEP (Steps 3-6)"
fi


# List the copied files for verification (only when not restoring from backup)
if [ "$STARTING_STEP" -le 2 ] || [ "$STARTING_STEP" -eq 7 ]; then
    log_message "Listing copied files:"
    ls -l 
    log_message "Processing sample: $SAMPLE"
else
    log_message "Skipping file listing - files restored from backup for starting step $STARTING_STEP"
fi

# Final verification of starting step before pipeline execution
log_message "Final verification - STARTING_STEP is set to: $STARTING_STEP"

# Function to verify STARTING_STEP hasn't changed
verify_starting_step() {
    local step_number="$1"
    log_message "Verifying STARTING_STEP before Step $step_number: STARTING_STEP=$STARTING_STEP"
    if [ "$STARTING_STEP" -ne "$step_number" ] && [ "$STARTING_STEP" -gt "$step_number" ]; then
        log_message "⚠️  Warning: STARTING_STEP ($STARTING_STEP) is greater than current step ($step_number). This step should be skipped."
    fi
}

# =============================================================================
# STEP 1: QUALITY CONTROL OF RAW READS
# =============================================================================
# FastQC provides comprehensive quality reports for raw sequencing data
# This helps identify potential issues before processing

echo ""
echo "┌─────────────────────────────────────────────────────────────────────────────────────┐"
echo "│                           🧬 STEP 1: QUALITY CONTROL OF RAW READS 🧬                │"
echo "│                           FastQC Analysis of Raw Sequencing Data                    │"
echo "└─────────────────────────────────────────────────────────────────────────────────────┘"
echo ""

log_message "Step 1 check: STARTING_STEP=$STARTING_STEP, condition: STARTING_STEP <= 1"
verify_starting_step 1
if [ "$STARTING_STEP" -le 1 ]; then
    start_step_timer "Step 1: FastQC on raw reads"
    if ! fastqc -t $NUM_THREADS *.fastq.gz ; then
        error_exit "FastQC failed"
    fi
    end_step_timer "Step 1: FastQC on raw reads"

    # Copy the fastqc results directly to RDM permanent storage
    rsync -rhivPt *.html "$RDM_QC_PRE/${SAMPLE}/" || error_exit "Failed to copy FastQC HTML results to RDM"
    rsync -rhivPt *.zip "$RDM_QC_PRE/${SAMPLE}/" || error_exit "Failed to copy FastQC ZIP results to RDM"

    # Remove the unneeded fastqc results from the $TMPDIR to save space
    rm *.html *.zip


    # No backup needed for FastQC on raw reads
else
    log_message "Step 1: FastQC on raw reads - SKIPPED (already completed)"
fi

# List files in TMPDIR for verification
log_message "Files in TMPDIR after Step 1:"
ls -l

# =============================================================================
# STEP 2: READ TRIMMING AND QUALITY FILTERING
# =============================================================================
# Trimmomatic removes adapters and low-quality bases
# This step is crucial for accurate alignment and variant calling

echo ""
echo "┌─────────────────────────────────────────────────────────────────────────────────────┐"
echo "│                           ✂️  STEP 2: READ TRIMMING AND QUALITY FILTERING ✂️           │"
echo "│                           Trimmomatic Adapter Removal and Quality Control             │"
echo "└─────────────────────────────────────────────────────────────────────────────────────┘"
echo ""

log_message "Step 2 check: STARTING_STEP=$STARTING_STEP, condition: STARTING_STEP <= 2"
verify_starting_step 2
if [ "$STARTING_STEP" -le 2 ]; then
    start_step_timer "Step 2: Read trimming with Trimmomatic"
    if ! java -Xmx${MEMORY} -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE \
        -threads $NUM_THREADS \
        *1.fastq.gz *2.fastq.gz \
        "${SAMPLE}_forward_paired.fastq.gz" "${SAMPLE}_forward_unpaired.fastq.gz" \
        "${SAMPLE}_reverse_paired.fastq.gz" "${SAMPLE}_reverse_unpaired.fastq.gz" \
        ILLUMINACLIP:$(basename "$ADAPTER_FILE"):${ILLUMINA_CLIP_SETTINGS} \
        LEADING:${TRIM_LEADING} TRAILING:${TRIM_TRAILING} \
        SLIDINGWINDOW:${SLIDING_WINDOW} AVGQUAL:${MIN_QUALITY} \
        MINLEN:${MIN_LENGTH}; then
        error_exit "Trimmomatic failed"
    fi
    end_step_timer "Step 2: Read trimming with Trimmomatic"

    # Trimmomatic parameters explained:
    # PE: Paired-end mode
    # ILLUMINACLIP: Remove Illumina adapters with specified settings
    # LEADING/TRAILING: Remove low-quality bases from read ends
    # SLIDINGWINDOW: Remove low-quality regions using sliding window
    # AVGQUAL: Minimum average quality score
    # MINLEN: Minimum read length after trimming

    # Removing the raw reads from the $TMPDIR to save space
    rm *1.fastq.gz *2.fastq.gz

    # Running FastQC on the trimmed reads to verify quality improvement
    start_step_timer "Step 2b: FastQC on trimmed reads"
    if ! fastqc -t $NUM_THREADS *.fastq.gz; then
        error_exit "FastQC on trimmed reads failed"
    fi
    end_step_timer "Step 2b: FastQC on trimmed reads"

    # Copy the fastqc results directly to RDM permanent storage
    rsync -rhivPt *.html "$RDM_QC_POST/${SAMPLE}/" || error_exit "Failed to copy FastQC HTML results to RDM"
    rsync -rhivPt *.zip "$RDM_QC_POST/${SAMPLE}/" || error_exit "Failed to copy FastQC ZIP results to RDM"

    # Remove the unneeded fastqc results from the $TMPDIR
    rm *.html *.zip

    # Create backup of Step 2 outputs (trimmed FASTQ files needed for Step 3)
    backup_current_state "02_read_trimming"
else
    log_message "Step 2: Read trimming - SKIPPED (already completed)"
fi

# List files in TMPDIR for verification
log_message "Files in TMPDIR after Step 2:"
ls -l

# =============================================================================
# STEP 3: READ ALIGNMENT TO REFERENCE GENOME
# =============================================================================
# BWA-MEM aligns trimmed reads to the reference genome
# This creates SAM files that are then converted to sorted BAM files

echo ""
echo "┌─────────────────────────────────────────────────────────────────────────────────────┐"
echo "│                           🎯 STEP 3: READ ALIGNMENT TO REFERENCE GENOME 🎯            │"
echo "│                           BWA-MEM Alignment and BAM File Generation                  │"
echo "└─────────────────────────────────────────────────────────────────────────────────────┘"
echo ""

log_message "Step 3 check: STARTING_STEP=$STARTING_STEP, condition: STARTING_STEP <= 3"
verify_starting_step 3
if [ "$STARTING_STEP" -le 3 ]; then
    start_step_timer "Step 3: BWA alignment"
    if ! bwa mem -M -t $NUM_THREADS -R "$READ_GROUP" \
        $REF_GENOME \
        "${SAMPLE}_forward_paired.fastq.gz" \
        "${SAMPLE}_reverse_paired.fastq.gz" \
        | samtools sort -o "${SAMPLE}_sorted.bam"; then
        error_exit "BWA alignment failed"
    fi
    end_step_timer "Step 3: BWA alignment"

    # BWA-MEM parameters explained:
    # -M: Mark shorter split hits as secondary (for Picard compatibility)
    # -t: Number of threads
    # -R: Read group information (essential for GATK)

    # Index the sorted BAM file for efficient access
    start_step_timer "Step 3b: BAM indexing"
    samtools index -@ $NUM_THREADS "${SAMPLE}_sorted.bam"
    end_step_timer "Step 3b: BAM indexing"

    # Remove all FASTQ files from TMPDIR to save space after alignment
    rm "${SAMPLE}_forward_paired.fastq.gz" "${SAMPLE}_reverse_paired.fastq.gz"
    rm "${SAMPLE}_forward_unpaired.fastq.gz" "${SAMPLE}_reverse_unpaired.fastq.gz"

    # Copy the mapped reads and the sorted BAM file directly to RDM permanent storage
    rsync -rhivPt "${SAMPLE}_sorted.bam" "$RDM_BAM/${SAMPLE}/" || error_exit "Failed to copy BAM file to RDM"
    rsync -rhivPt "${SAMPLE}_sorted.bam.bai" "$RDM_BAM/${SAMPLE}/" || error_exit "Failed to copy BAM index to RDM"

    log_message "BWA alignment and indexing completed successfully"

    # Create backup of Step 3 outputs
    backup_current_state "03_bwa_alignment"
else
    log_message "Step 3: BWA alignment - SKIPPED (already completed)"
fi



# List files in TMPDIR for verification
log_message "Files in TMPDIR after Step 3:"
ls -l
echo "-----------------------------"
echo "."
echo "."
echo "."
echo "."
echo "."
echo "."
echo "."
echo "."

# =============================================================================
# STEP 4: DUPLICATE MARKING
# =============================================================================
# MarkDuplicates identifies and marks PCR duplicates
# These are reads that appear to be identical due to PCR amplification
# Marking them helps prevent false variant calls

echo ""
echo "┌─────────────────────────────────────────────────────────────────────────────────────┐"
echo "│                           🔍 STEP 4: DUPLICATE MARKING 🔍                            │"
echo "│                           GATK MarkDuplicates for PCR Duplicate Removal              │"
echo "└─────────────────────────────────────────────────────────────────────────────────────┘"
echo ""

log_message "Step 4 check: STARTING_STEP=$STARTING_STEP, condition: STARTING_STEP <= 4"
verify_starting_step 4
if [ "$STARTING_STEP" -le 4 ]; then
    start_step_timer "Step 4: Duplicate marking"
    if ! gatk --java-options "-Xmx${MEMORY}" MarkDuplicates \
        -I "${SAMPLE}_sorted.bam" \
        -O "${SAMPLE}_dedup.bam" \
        -M "${SAMPLE}_dedup_metrics_file.txt"; then
        error_exit "MarkDuplicates failed"
    fi
    end_step_timer "Step 4: Duplicate marking"

    # Index the deduplicated BAM file
    start_step_timer "Step 4b: Deduplicated BAM indexing"
    samtools index -@ $NUM_THREADS "${SAMPLE}_dedup.bam"
    end_step_timer "Step 4b: Deduplicated BAM indexing"

    # Remove the unneeded files from the $TMPDIR
    rm "${SAMPLE}_sorted.bam" "${SAMPLE}_sorted.bam.bai"

    # Copy deduplicated BAM and index directly to RDM permanent storage
    rsync -rhivPt "${SAMPLE}_dedup.bam" "${SAMPLE}_dedup.bam.bai" "$RDM_BAM/${SAMPLE}/" || error_exit "Failed to sync deduplicated BAM to RDM"

    # Sync MarkDuplicates metrics to RDM Metrics directory
    rsync -rhivPt "${SAMPLE}_dedup_metrics_file.txt" "$RDM_METRICS/${SAMPLE}/" || error_exit "Failed to sync metrics file to RDM Metrics directory"

    # Remove obsolete sorted BAM files from backup before creating new backup snapshot
    prune_backup_files "post-markduplicates" "${SAMPLE}_sorted.bam" "${SAMPLE}_sorted.bam.bai"

    # Create backup of Step 4 outputs
    backup_current_state "04_mark_duplicates"
else
    log_message "Step 4: Duplicate marking - SKIPPED (already completed)"
fi


# List files in TMPDIR for verification
log_message "Files in TMPDIR after Step 4:"
ls -l

# =============================================================================
# STEP 5: BASE QUALITY SCORE RECALIBRATION (BQSR)
# =============================================================================
# BQSR corrects systematic errors in base quality scores
# This improves the accuracy of variant calling
# Uses known variant sites to train the recalibration model

echo ""
echo "┌─────────────────────────────────────────────────────────────────────────────────────┐"
echo "│                           📊 STEP 5: BASE QUALITY SCORE RECALIBRATION (BQSR) 📊      │"
echo "│                           GATK BaseRecalibrator and ApplyBQSR                       │"
echo "└─────────────────────────────────────────────────────────────────────────────────────┘"
echo ""

log_message "Step 5 check: STARTING_STEP=$STARTING_STEP, condition: STARTING_STEP <= 5"
verify_starting_step 5
if [ "$STARTING_STEP" -le 5 ]; then
    # Build the recalibration table
    start_step_timer "Step 5a: Base recalibration table generation"
    if ! gatk --java-options "-Xmx${MEMORY}" BaseRecalibrator \
        -R $REF_GENOME \
        -I "${SAMPLE}_dedup.bam" \
        --known-sites "$KNOWN_SITES_PATH" \
        -O "${SAMPLE}_recal.table"; then
        error_exit "BaseRecalibrator failed"
    fi
    end_step_timer "Step 5a: Base recalibration table generation"

    # Apply the recalibration table to create final BAM
    start_step_timer "Step 5b: Apply base recalibration"
    if ! gatk --java-options "-Xmx${MEMORY}" ApplyBQSR \
        -R $REF_GENOME \
        -I "${SAMPLE}_dedup.bam" \
        --bqsr-recal-file "${SAMPLE}_recal.table" \
        -O "${SAMPLE}_recal.bam"; then
        error_exit "ApplyBQSR failed"
    fi
    end_step_timer "Step 5b: Apply base recalibration"

    # Index the recalibrated BAM file
    start_step_timer "Step 5c: Recalibrated BAM indexing"
    samtools index -@ $NUM_THREADS "${SAMPLE}_recal.bam"
    end_step_timer "Step 5c: Recalibrated BAM indexing"

    log_message "Base recalibration completed successfully"

    # Copy the recalibrated BAM file and index directly to RDM permanent storage
    rsync -rhivPt "${SAMPLE}_recal.bam" "$RDM_BAM/${SAMPLE}/" || error_exit "Failed to copy recalibrated BAM to RDM"
    rsync -rhivPt "${SAMPLE}_recal.bam.bai" "$RDM_BAM/${SAMPLE}/" || error_exit "Failed to copy BAM index to RDM"

    # Remove obsolete deduplicated BAM from backup before creating new backup snapshot
    prune_backup_files "post-bqsr" "${SAMPLE}_dedup.bam" "${SAMPLE}_dedup.bam.bai"

    # Create backup of Step 5 outputs
    backup_current_state "05_base_recalibration"

    # Remove the unneeded files from the $TMPDIR
    rm "${SAMPLE}_dedup.bam" "${SAMPLE}_dedup.bam.bai"
else
    log_message "Step 5: Base quality score recalibration - SKIPPED (already completed)"
fi

# List files in TMPDIR for verification
log_message "Files in TMPDIR after Step 5:"
ls -l

# =============================================================================
# STEP 6: VARIANT CALLING WITH HAPLOTYPE CALLER
# =============================================================================
# HaplotypeCaller performs local de novo assembly of haplotypes
# This is the most accurate variant caller in GATK
# GVCF output allows for joint genotyping in the next stage

echo ""
echo "┌─────────────────────────────────────────────────────────────────────────────────────┐"
echo "│                           🧬 STEP 6: VARIANT CALLING WITH HAPLOTYPE CALLER 🧬      │"
echo "│                           GATK HaplotypeCaller for GVCF Generation                 │"
echo "└─────────────────────────────────────────────────────────────────────────────────────┘"
echo ""

log_message "Step 6 check: STARTING_STEP=$STARTING_STEP, condition: STARTING_STEP <= 6"
verify_starting_step 6
if [ "$STARTING_STEP" -le 6 ]; then
    start_step_timer "Step 6: Variant calling with HaplotypeCaller"
    if ! gatk --java-options "-Xmx${MEMORY}" HaplotypeCaller \
        -R $REF_GENOME \
        -I "${SAMPLE}_recal.bam" \
        -O "${SAMPLE}_raw.g.vcf.gz" \
        -ERC GVCF \
        --native-pair-hmm-threads $NUM_THREADS; then
        error_exit "HaplotypeCaller failed"
    fi
    end_step_timer "Step 6: Variant calling with HaplotypeCaller"

    # HaplotypeCaller parameters explained:
    # --java-options "-Xmx24g": Allocate 24GB memory to Java
    # -R: Reference genome
    # -I: Input BAM file
    # -O: Output GVCF file
    # -ERC GVCF: Emit reference confidence scores (required for joint genotyping)
    # --native-pair-hmm-threads: Number of threads for PairHMM algorithm

    # Sync the raw GVCF file to RDM permanent storage immediately after HaplotypeCaller
    # This protects the expensive HaplotypeCaller output in case of system failure
    start_step_timer "Step 6b: GVCF sync to RDM"
    rsync -rhivPt "${SAMPLE}_raw.g.vcf.gz" "$RDM_VCF/" || error_exit "Failed to sync GVCF file to RDM"
    end_step_timer "Step 6b: GVCF sync to RDM"

    # Index the GVCF file for efficient access and GenomicsDBImport compatibility
    start_step_timer "Step 6c: GVCF indexing"
    if ! gatk IndexFeatureFile --input "${SAMPLE}_raw.g.vcf.gz"; then
        error_exit "GVCF indexing failed"
    fi
    end_step_timer "Step 6c: GVCF indexing"

    # Sync the index file to RDM permanent storage
    start_step_timer "Step 6d: GVCF index sync to RDM"
    rsync -rhivPt "${SAMPLE}_raw.g.vcf.gz.tbi" "$RDM_VCF/" || error_exit "Failed to sync GVCF index to RDM"
    end_step_timer "Step 6d: GVCF index sync to RDM"

    # Create backup of Step 6 outputs (comprehensive backup)
    backup_current_state "06_variant_calling"
else
    log_message "Step 6: Variant calling - SKIPPED (already completed)"
fi

# =============================================================================
# STEP 7: GENOTYPE GVCF TO VCF CONVERSION
# =============================================================================
# GenotypeGVCFs converts the GVCF file to a standard VCF with actual genotype calls
# This step is separate from Step 6 to avoid rerunning expensive HaplotypeCaller
# if only the genotyping step fails

echo ""
echo "┌─────────────────────────────────────────────────────────────────────────────────────┐"
echo "│                           🧬 STEP 7: GENOTYPE GVCF TO VCF CONVERSION 🧬             │"
echo "│                           GATK GenotypeGVCFs for VCF Generation                     │"
echo "└─────────────────────────────────────────────────────────────────────────────────────┘"
echo ""

log_message "Step 7 check: STARTING_STEP=$STARTING_STEP, condition: STARTING_STEP <= 7"
verify_starting_step 7
if [ "$STARTING_STEP" -le 7 ]; then
    start_step_timer "Step 7: GenotypeGVCFs"
    if ! gatk --java-options "-Xmx${MEMORY}" GenotypeGVCFs \
        -R $REF_GENOME \
        -V "${SAMPLE}_raw.g.vcf.gz" \
        -O "${SAMPLE}_genotyped.vcf.gz"; then
        error_exit "GenotypeGVCFs failed"
    fi
    end_step_timer "Step 7: GenotypeGVCFs"

    # GenotypeGVCFs parameters explained:
    # --java-options "-Xmx24g": Allocate 24GB memory to Java
    # -R: Reference genome
    # -V: Input GVCF file
    # -O: Output VCF file with genotype calls

    # Copy the genotyped VCF file directly to RDM permanent storage
    rsync -rhivPt "${SAMPLE}_genotyped.vcf.gz" "$RDM_VCF/" || error_exit "Failed to copy genotyped VCF file to RDM"

else
    log_message "Step 7: GenotypeGVCFs - SKIPPED (already completed)"
fi

# List files in TMPDIR for verification
log_message "Files in TMPDIR after Step 7:"
ls -l

# =============================================================================
# PIPELINE COMPLETION WITH TIMING SUMMARY
# =============================================================================
# At this point, the per-sample processing is complete
# Both GVCF and VCF files are ready for analysis

log_message "Pipeline Stage 1A completed successfully for sample: $SAMPLE"
log_message "GVCF file saved to: $RDM_VCF/${SAMPLE}_raw.g.vcf.gz"
log_message "Genotyped VCF file saved to: $RDM_VCF/${SAMPLE}_genotyped.vcf.gz"
log_message "Analysis-ready BAM file saved to: $RDM_BAM/${SAMPLE}/${SAMPLE}_recal.bam"
log_message "All backups saved to: $BACKUP_DIR"
log_message "RDM permanent storage updated with all outputs"

# Display comprehensive timing summary
show_total_time

# =============================================================================
# FINAL CLEANUP: REMOVE ALL BACKUP FILES
# =============================================================================
# After successful pipeline completion, remove all backup files to save storage space
# All outputs are safely stored in RDM permanent storage

log_message "Pipeline completed successfully - removing backup files to save storage space"
if [ -d "$BACKUP_DIR" ]; then
    log_message "Removing backup directory: $BACKUP_DIR"
    rm -rf "$BACKUP_DIR"
    log_message "Backup cleanup completed"
else
    log_message "No backup directory found to clean up"
fi




