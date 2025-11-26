#!/bin/bash -l
#SBATCH --job-name=Apple_GATK_1B_KH
#SBATCH --ntasks=1
#SBATCH --account=a_qaafi_cas
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --cpus-per-task=6
#SBATCH --time=200:00:00
#SBATCH --mem=36G
#SBATCH --array=0-17
#SBATCH -o /scratch/user/uqpha1/logs/1B_CombineGVCFs/%A_%a_%x_%j.output
#SBATCH -e /scratch/user/uqpha1/logs/1B_CombineGVCFs/%A_%a_%x_%j.error

# =============================================================================
# APPLE GATK PIPELINE - STAGE 1B: COMBINE INDIVIDUAL GVCFs
# =============================================================================
# AUTHORS: Phu Khang Ha
# DATE: 1/10/2025
#
# PIPELINE OVERVIEW:
# This script combines individual GVCF files from Stage 1A into consolidated
# VCF files using GenomicsDBImport and GenotypeGVCFs.
# Since this step is very intensive, the script will utilise /scratch as storage for the intermediate files.
# 
# WORKFLOW STEPS:
# 1. Extract RDM folder path and create scratch working directory
# 2. Create tar.gz archive of all GVCF files from RDM/5.Individual_VCF (Task 0 only)
# 3. Extract GVCF files to scratch working directory (Task 0 only)
# 4. Run GenomicsDBImport for assigned chromosome (each array task)
# 5. Run GenotypeGVCFs for assigned chromosome (each array task)
# 6. Copy consolidated VCF files back to RDM/6.Consolidated_GVCF (each array task)
#
# INPUT: RDM base path containing individual GVCF files
# OUTPUT: Consolidated VCF files per chromosome
# =============================================================================

# Get command line arguments
rdm_base_path=$1         # RDM base path containing the standard folder structure

# Extract dataset name from RDM path for scratch directory naming
today_dataset_folder=$(basename "$rdm_base_path")

# Define main working directory
readonly workdir="/scratch/user/uqpha1/${today_dataset_folder}_1B"

# Validate command line arguments
if [ -z "$rdm_base_path" ]; then
    echo "Error: Missing required arguments"
    echo "Usage: sbatch script.sh <rdm_base_path>"
    echo "Example: sbatch script.sh /QRISdata/Q8367/NCBI"
    echo ""
    echo "RDM base path should contain the standard folder structure:"
    echo "  /QRISdata/Q8367/NCBI/"
    echo "  ├── 5.Individual_VCF/     # Individual GVCF files from Stage 1A"
    echo "  └── 6.Consolidated_GVCF/  # Output directory for consolidated VCFs"
    exit 1
fi

# Validate RDM base path exists and is a directory
if [ ! -d "$rdm_base_path" ]; then
    echo "Error: RDM base path '$rdm_base_path' does not exist or is not a directory"
    exit 1
fi

# Validate that RDM base path contains the expected folder structure
if [ ! -d "$rdm_base_path/5.Individual_VCF" ]; then
    echo "Error: RDM base path '$rdm_base_path' does not contain expected folder structure"
    echo "Expected: $rdm_base_path/5.Individual_VCF"
    echo "Please ensure the RDM path contains individual GVCF files from Stage 1A"
    exit 1
fi

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

# Logging function - adds timestamps to all messages
log_message() {
    echo ""
    echo "=================================================================================="
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
    echo ""
}

# Error handling function - exits with error message
error_exit() {
    local error_message="$1"
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    
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
    echo "    🕐 Error Time: $timestamp"
    echo "    🔍 Error Message: $error_message"
    echo ""
    echo "████████████████████████████████████████████████████████████████████████████████████████"
    echo "████████████████████████████████████████████████████████████████████████████████████████"
    echo "████████████████████████████████████████████████████████████████████████████████████████"
    echo ""
    echo ""
    
    echo "Error: $error_message" >&2
    exit 1
}

# =============================================================================
# MODULE LOADING SECTION
# =============================================================================

module load gatk/4.3.0.0-gcccore-11.3.0-java-11

# =============================================================================
# PIPELINE PARAMETERS AND CONSTANTS
# =============================================================================

readonly NUM_THREADS=12                         # Number of CPU threads for parallel processing
readonly MEMORY="36G"                          # Memory allocation for Java-based tools
readonly JAVA_TMP_SIZE="4G"                    # Java temporary directory size

# Apple genome chromosomes (18 chromosomes: Chr00-Chr17)
readonly CHROMOSOMES=(00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17)

# Get current chromosome from array task ID
CURRENT_CHR="${CHROMOSOMES[$SLURM_ARRAY_TASK_ID]}"
CURRENT_CONTIG="Chr${CURRENT_CHR}"

# =============================================================================
# REFERENCE FILES AND PATHS
# =============================================================================

readonly REF_GENOME_DIR="/QRISdata/Q8367/Reference_Genome"   # Directory containing reference genome
readonly REF_GENOME="GDDH13_1-1_formatted.fasta"            # Apple reference genome (GDDH13 v1.1)

# RDM (Permanent) storage paths
readonly RDM_BASE="$rdm_base_path"
readonly RDM_INDIVIDUAL_VCF="${RDM_BASE}/5.Individual_VCF"
readonly RDM_CONSOLIDATED="${RDM_BASE}/6.Consolidated_GVCF"

# Shared reference directory in scratch (similar to 1A approach)
readonly SHARED_REF_GENOME_DIR="${workdir}/9.Reference_genome_dir"

# =============================================================================
# DIRECTORY CREATION AND INITIALIZATION
# =============================================================================

log_message "Initializing Apple GATK Pipeline Stage 1B: Combine GVCFs"
log_message "RDM base path: $rdm_base_path"
log_message "Working directory: $workdir"
log_message "Array task ID: $SLURM_ARRAY_TASK_ID"
log_message "Processing chromosome: $CURRENT_CONTIG"

# Create working directory (only by task 0)
if [ "$SLURM_ARRAY_TASK_ID" -eq 0 ]; then
    mkdir -p "$workdir" || error_exit "Failed to create working directory: $workdir"
    mkdir -p "/scratch/user/uqpha1/logs/1B_CombineGVCFs" || error_exit "Failed to create log directory"
    mkdir -p "$RDM_CONSOLIDATED" || error_exit "Failed to create consolidated VCF directory: $RDM_CONSOLIDATED"
    mkdir -p "$SHARED_REF_GENOME_DIR" || error_exit "Failed to create shared reference directory: $SHARED_REF_GENOME_DIR"
fi

# Wait for task 0 to create directories
if [ "$SLURM_ARRAY_TASK_ID" -gt 0 ]; then
    sleep 10
fi

# Change to working directory
cd "$workdir" || error_exit "Failed to change to working directory: $workdir"

echo ""
echo ""
echo "████████████████████████████████████████████████████████████████████████████████████████"
echo "████████████████████████████████████████████████████████████████████████████████████████"
echo "████████████████████████████████████████████████████████████████████████████████████████"
echo ""
echo "                                    🚀 APPLE GATK PIPELINE 1B STARTING 🚀"
echo ""
echo "████████████████████████████████████████████████████████████████████████████████████████"
echo "████████████████████████████████████████████████████████████████████████████████████████"
echo "████████████████████████████████████████████████████████████████████████████████████████"
echo ""
echo "    📁 RDM Base Path: $rdm_base_path"
echo "    🏠 Working Directory: $workdir"
echo "    🧬 Processing Chromosome: $CURRENT_CONTIG"
echo "    🔢 Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "    🕐 Started: $(date)"
echo ""
echo "████████████████████████████████████████████████████████████████████████████████████████"
echo "████████████████████████████████████████████████████████████████████████████████████████"
echo "████████████████████████████████████████████████████████████████████████████████████████"
echo ""
echo ""

# =============================================================================
# TIMING FUNCTIONS
# =============================================================================

# Global timing variables
PIPELINE_START_TIME=$(date +%s)
STEP_START_TIME=0
STEP_END_TIME=0
CURRENT_STEP_NAME=""
STEP_TIMES_FILE="${workdir}/step_times.log"

# Function to log step timing to file
log_step_time() {
    local step_name="$1"
    local timestamp="$2"
    local status="$3"
    echo "[$timestamp] $status: $step_name" >> "$STEP_TIMES_FILE"
}

# Function to start timing a step
start_step_timer() {
    local step_name="$1"
    STEP_START_TIME=$(date +%s)
    CURRENT_STEP_NAME="$step_name"
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    
    echo ""
    echo "╔══════════════════════════════════════════════════════════════════════════════════╗"
    echo "║                           STARTING: $step_name"
    echo "║                           Started at: $timestamp"
    echo "╚══════════════════════════════════════════════════════════════════════════════════╝"
    echo ""
    
    # Log step start to file
    log_step_time "$step_name" "$timestamp" "STARTED"
}

# Function to end timing a step
end_step_timer() {
    local step_name="$1"
    STEP_END_TIME=$(date +%s)
    local duration=$((STEP_END_TIME - STEP_START_TIME))
    local hours=$((duration / 3600))
    local minutes=$(((duration % 3600) / 60))
    local seconds=$((duration % 60))
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    
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

# Initialize timing log file
log_message "Initializing timing log: $STEP_TIMES_FILE"
log_step_time "PIPELINE_START" "$(date '+%Y-%m-%d %H:%M:%S')" "PIPELINE_INITIALIZED"

# =============================================================================
# STEP 1: CREATE TAR.GZ ARCHIVE OF GVCF FILES (Task 0 only)
# =============================================================================

if [ "$SLURM_ARRAY_TASK_ID" -eq 0 ]; then
    start_step_timer "Step 1: Create tar.gz archive of GVCF files"

    # Check if GVCF files exist
    gvcf_count=$(find "$RDM_INDIVIDUAL_VCF" -name "*_raw.g.vcf.gz" | wc -l)
    tbi_count=$(find "$RDM_INDIVIDUAL_VCF" -name "*_raw.g.vcf.gz.tbi" | wc -l)

    if [ "$gvcf_count" -eq 0 ]; then
        error_exit "No GVCF files found in $RDM_INDIVIDUAL_VCF"
    fi

    log_message "Found $gvcf_count GVCF files and $tbi_count TBI index files to process"

    # Create tar.gz archive with compression
    archive_name="gvcf_files_${today_dataset_folder}_$(date +%Y%m%d_%H%M%S).tar.gz"
    log_message "Creating compressed archive: $archive_name"

    # Create archive with gzip compression level 9 (maximum compression)
    if ! tar -czf "$archive_name" -C "$RDM_INDIVIDUAL_VCF" .; then
        error_exit "Failed to create tar.gz archive of GVCF files"
    fi

    # Verify archive was created
    if [ ! -f "$archive_name" ]; then
        error_exit "Archive file was not created: $archive_name"
    fi

    archive_size=$(du -h "$archive_name" | cut -f1)
    log_message "Archive created successfully: $archive_name (Size: $archive_size)"

    end_step_timer "Step 1: Create tar.gz archive of GVCF files"
else
    log_message "Step 1: Skipped (Task 0 only) - waiting for archive creation"
fi

# =============================================================================
# STEP 2: EXTRACT GVCF FILES TO WORKING DIRECTORY (Task 0 only)
# =============================================================================

if [ "$SLURM_ARRAY_TASK_ID" -eq 0 ]; then
    start_step_timer "Step 2: Extract GVCF files to working directory"

    log_message "Extracting GVCF files from archive"

    if ! tar -xzf "$archive_name"; then
        error_exit "Failed to extract GVCF files from archive"
    fi

    # Verify extraction
    extracted_gvcf_count=$(find . -maxdepth 1 -name "*_raw.g.vcf.gz" | wc -l)
    extracted_tbi_count=$(find . -maxdepth 1 -name "*_raw.g.vcf.gz.tbi" | wc -l)

    if [ "$extracted_gvcf_count" -ne "$gvcf_count" ]; then
        error_exit "Extraction failed: expected $gvcf_count GVCF files, found $extracted_gvcf_count"
    fi

    if [ "$extracted_tbi_count" -ne "$tbi_count" ]; then
        error_exit "Extraction failed: expected $tbi_count TBI files, found $extracted_tbi_count"
    fi

    log_message "Successfully extracted $extracted_gvcf_count GVCF files and $extracted_tbi_count TBI index files"

    # Remove archive to save space
    rm "$archive_name"
    log_message "Removed archive file to save space"

    end_step_timer "Step 2: Extract GVCF files to working directory"
else
    log_message "Step 2: Skipped (Task 0 only) - waiting for file extraction"
fi

# =============================================================================
# STEP 3: SETUP SHARED REFERENCE GENOME (Task 0 only)
# =============================================================================

if [ "$SLURM_ARRAY_TASK_ID" -eq 0 ]; then
    start_step_timer "Step 3: Setup shared reference genome"

    log_message "Setting up shared reference genome directory"

    # Check if reference genome already exists in shared directory
    if [ -f "${SHARED_REF_GENOME_DIR}/${REF_GENOME}" ]; then
        log_message "Reference genome already exists in shared directory: ${SHARED_REF_GENOME_DIR}"
    else
        log_message "Copying reference genome from RDM to shared directory"
        
        if ! rsync -rhivPt "${REF_GENOME_DIR}/" "${SHARED_REF_GENOME_DIR}/"; then
            error_exit "Failed to copy reference genome files to shared directory"
        fi

        # Verify reference genome was copied
        if [ ! -f "${SHARED_REF_GENOME_DIR}/${REF_GENOME}" ]; then
            error_exit "Reference genome file not found: ${SHARED_REF_GENOME_DIR}/${REF_GENOME}"
        fi

        log_message "Reference genome copied successfully to shared directory"
    fi

    end_step_timer "Step 3: Setup shared reference genome"
else
    log_message "Step 3: Skipped (Task 0 only) - waiting for shared reference genome setup"
fi

# =============================================================================
# STEP 4: CREATE MISSING INDEX FILES (Task 0 only)
# =============================================================================

if [ "$SLURM_ARRAY_TASK_ID" -eq 0 ]; then
    start_step_timer "Step 4: Create missing index files"

    log_message "Checking for missing .tbi index files"

    missing_indexes=0
    for gvcf_file in ./*_raw.g.vcf.gz; do
        if [ -f "$gvcf_file" ]; then
            tbi_file="${gvcf_file}.tbi"
            if [ ! -f "$tbi_file" ]; then
                log_message "Creating missing index: $tbi_file"
                if ! gatk IndexFeatureFile --input "$gvcf_file"; then
                    error_exit "Failed to create index: $tbi_file"
                fi
                ((missing_indexes++))
            fi
        fi
    done

    if [ $missing_indexes -eq 0 ]; then
        log_message "All GVCF files already have .tbi indexes (extracted from archive)"
    else
        log_message "Created $missing_indexes missing .tbi files"
    fi

    end_step_timer "Step 4: Create missing index files"
else
    log_message "Step 4: Skipped (Task 0 only) - waiting for index creation"
fi

# =============================================================================
# STEP 5: RUN GENOMICSDBIMPORT FOR ASSIGNED CHROMOSOME
# =============================================================================

start_step_timer "Step 5: Run GenomicsDBImport for $CURRENT_CONTIG"

# Create GenomicsDB output directory
mkdir -p "genomicsdb_output"

# Create symlinks to shared reference genome files
log_message "Creating symlinks to shared reference genome files"
ln -sf "${SHARED_REF_GENOME_DIR}"/* .

# Set reference genome path to shared directory
readonly REF_GENOME_PATH="${SHARED_REF_GENOME_DIR}/${REF_GENOME}"

# Get list of GVCF files
gvcf_files=$(find . -maxdepth 1 -name "*_raw.g.vcf.gz" | sort | tr '\n' ' ')
log_message "GVCF files to process: $gvcf_files"

log_message "Processing $CURRENT_CONTIG with GenomicsDBImport..."

chr_output_dir="genomicsdb_output/${CURRENT_CONTIG}_gdb"
chr_log="genomicsdb_${CURRENT_CONTIG}.log"

# Remove existing workspace if it exists
if [ -d "$chr_output_dir" ]; then
    log_message "Removing existing workspace: $chr_output_dir"
    rm -rf "$chr_output_dir"
fi

# Run GenomicsDBImport for this chromosome
if ! gatk --java-options "-Djava.io.tmpdir=$TMPDIR -Xms${JAVA_TMP_SIZE} -Xmx${MEMORY} -XX:ParallelGCThreads=4" \
    GenomicsDBImport \
    --verbosity INFO \
    --genomicsdb-workspace-path "$chr_output_dir" \
    --overwrite-existing-genomicsdb-workspace true \
    --reference "$REF_GENOME_PATH" \
    $(for file in $gvcf_files; do echo "--variant $file"; done) \
    --tmp-dir "$TMPDIR" \
    --max-num-intervals-to-import-in-parallel 6 \
    --intervals "$CURRENT_CONTIG" \
    --batch-size 100 \
    --reader-threads 4 \
    --genomicsdb-shared-posixfs-optimizations true \
    --consolidate true \
    --merge-input-intervals true \
    2>&1 | tee "$chr_log"; then
    
    error_exit "Failed to process $CURRENT_CONTIG with GenomicsDBImport"
else
    log_message "Successfully processed $CURRENT_CONTIG with GenomicsDBImport"
fi

end_step_timer "Step 5: Run GenomicsDBImport for $CURRENT_CONTIG"

# =============================================================================
# STEP 6: RUN GENOTYPEGVCFS FOR ASSIGNED CHROMOSOME
# =============================================================================

start_step_timer "Step 6: Run GenotypeGVCFs for $CURRENT_CONTIG"

log_message "Processing $CURRENT_CONTIG with GenotypeGVCFs..."

genomicsdb_path="genomicsdb_output/${CURRENT_CONTIG}_gdb"
output_vcf="${CURRENT_CONTIG}.vcf.gz"
chr_log="genotype_${CURRENT_CONTIG}.log"

# Check if GenomicsDB exists for this chromosome
if [ ! -d "$genomicsdb_path" ]; then
    error_exit "GenomicsDB not found for $CURRENT_CONTIG: $genomicsdb_path"
fi

# Check if output VCF already exists
if [ -f "$output_vcf" ]; then
    log_message "$CURRENT_CONTIG already processed: $output_vcf exists, skipping"
else
    # Run GenotypeGVCFs for this chromosome
    if ! gatk --java-options "-Djava.io.tmpdir=$TMPDIR -Xms${JAVA_TMP_SIZE} -Xmx${MEMORY} -XX:ParallelGCThreads=4" \
        GenotypeGVCFs \
        -R "$REF_GENOME_PATH" \
        -V "gendb://${genomicsdb_path}" \
        -O "$output_vcf" \
        2>&1 | tee "$chr_log"; then
        
        error_exit "Failed to process $CURRENT_CONTIG with GenotypeGVCFs"
    else
        log_message "Successfully processed $CURRENT_CONTIG with GenotypeGVCFs"
    fi
fi

end_step_timer "Step 6: Run GenotypeGVCFs for $CURRENT_CONTIG"

# =============================================================================
# STEP 7: COPY CONSOLIDATED VCF FILE TO RDM
# =============================================================================

start_step_timer "Step 7: Copy consolidated VCF file to RDM"

log_message "Copying consolidated VCF file to RDM: $RDM_CONSOLIDATED"

vcf_file="${CURRENT_CONTIG}.vcf.gz"

if [ -f "$vcf_file" ]; then
    if ! rsync -rhivPt "$vcf_file" "$RDM_CONSOLIDATED/"; then
        error_exit "Failed to copy VCF file: $vcf_file"
    fi
    log_message "Successfully copied: $vcf_file"
else
    error_exit "VCF file not found: $vcf_file"
fi

end_step_timer "Step 7: Copy consolidated VCF file to RDM"

# =============================================================================
# PIPELINE COMPLETION
# =============================================================================

# Calculate total pipeline time
pipeline_end_time=$(date +%s)
total_duration=$((pipeline_end_time - PIPELINE_START_TIME))
hours=$((total_duration / 3600))
minutes=$(((total_duration % 3600) / 60))
seconds=$((total_duration % 60))

log_message "Pipeline Stage 1B completed successfully for $CURRENT_CONTIG"
log_message "Consolidated VCF file saved to: $RDM_CONSOLIDATED/${CURRENT_CONTIG}.vcf.gz"
log_message "Working directory: $workdir"

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
echo "    📊 Dataset: $today_dataset_folder"
echo "    📁 RDM Base Path: $rdm_base_path"
echo "    🏠 Working Directory: $workdir"
echo "    🧬 Processed Chromosome: $CURRENT_CONTIG"
echo "    🔢 Array Task ID: $SLURM_ARRAY_TASK_ID"
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
echo "    📁 Consolidated VCF Location: $RDM_CONSOLIDATED"
echo "    📁 Step Times Log: $STEP_TIMES_FILE"
echo ""
echo "████████████████████████████████████████████████████████████████████████████████████████"
echo "████████████████████████████████████████████████████████████████████████████████████████"
echo "████████████████████████████████████████████████████████████████████████████████████████"
echo ""
echo ""

# Log final timing
log_step_time "PIPELINE_COMPLETE" "$(date '+%Y-%m-%d %H:%M:%S')" "EXIT_CODE_0"

log_message "Apple GATK Pipeline Stage 1B completed successfully!"
