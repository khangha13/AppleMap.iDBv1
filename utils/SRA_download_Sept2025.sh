#!/bin/bash -l
#SBATCH --job-name=K_download
#SBATCH --ntasks=1
#SBATCH --account=a_qaafi_cas
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --time=150:00:00
#SBATCH --mem=6G
#SBATCH --array=0-1
#SBATCH -o /scratch/user/uqpha1/logs/%x_%j_%a.output
#SBATCH -e /scratch/user/uqpha1/logs/%x_%j_%a.error

# =============================================================================
# SRA DOWNLOAD PIPELINE - ARRAY JOB VERSION
# =============================================================================
# Script to download SRA samples from NCBI using prefetch, validate with vdb-validate, 
# and convert to FASTQ with fasterq-dump.
# Each SLURM array task processes one sample from the sample list.
# 
# The sample_list path is hardcoded in the script.
# Update the sample_list path variable as needed for your run.
# 
# The sample_list should contain one SRA accession per line.
# =============================================================================

start_ts=$(date +%s)
# Set the sample list path (update this path as needed)
# Check for required argument: sample list file
if [ $# -ne 1 ]; then
    echo "Usage: $0 <sample_list.txt>"
    exit 1
fi

sample_list="$1"


# Get the sample ID from the sample list based on array task ID
echo "SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID"
SAMPLE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$sample_list")

if [ -z "$SAMPLE" ]; then
    echo "Error: No sample found for task ID ${SLURM_ARRAY_TASK_ID} (line $((SLURM_ARRAY_TASK_ID + 1)))"
    exit 1
fi

# =============================================================================
# DIRECTORY SETUP
# =============================================================================
# Set up working directory
DATE=$(date +%Y%m%d)
scratch_dir="/scratch/user/uqpha1/${DATE}_NCBI_sample"
if [ ! -d "$scratch_dir" ]; then
    mkdir -p "$scratch_dir"
    echo "Created directory: $scratch_dir"
else
    echo "Directory already exists: $scratch_dir"
fi

echo "Processing sample: $SAMPLE"
echo "Downloading $SAMPLE to $scratch_dir"
cd "$scratch_dir"

# =============================================================================
# DOWNLOAD AND PROCESSING
# =============================================================================

#Load module sra-toolkit
module purge
module load sra-toolkit

# Download SRA file with prefetch
echo "Downloading $SAMPLE with prefetch..."
prefetch --max-size 200GB -v -p --output-directory "$scratch_dir" "$SAMPLE"

# Validate SRA file
echo "Validating $SAMPLE"
vdb-validate "$scratch_dir/$SAMPLE"

# Convert to FASTQ using fasterq-dump, using $TMPDIR for temp files
echo "Running fasterq-dump for sample $SAMPLE"
fasterq-dump "$SAMPLE" -m 5GB -v -p -t "$TMPDIR"
# Compress FASTQ files
echo "Compressing fastq files for sample $SAMPLE"
gzip "$scratch_dir/$SAMPLE"*.fastq
# Test the integrity of the compressed FASTQ files
echo "Testing gzip integrity for sample $SAMPLE"
if ! gzip -t -v "$scratch_dir/$SAMPLE"*.fastq.gz; then
    echo "Error: gzip integrity test failed for sample $SAMPLE"
    exit 1
else
    echo "gzip integrity test passed for sample $SAMPLE"
fi


echo "Copying compressed files to final location"
rsync -rhivPt "$scratch_dir/$SAMPLE"*.fastq.gz /QRISdata/Q8367/WGS_Reference_Panel/NCBI/1.FASTQ/.

echo "Download and conversion complete for $SAMPLE"
echo "FASTQ files are in $scratch_dir"

end_ts=$(date +%s)
echo "Runtime (s): $((end_ts - start_ts))"