#!/bin/bash -l
#SBATCH --job-name=bam_2x_downsample
#SBATCH --ntasks=1
#SBATCH --account=a_qaafi_cas
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00
#SBATCH --mem=8G
#SBATCH -o /scratch/user/uqpha1/logs/%x_%j_%a.output
#SBATCH -e /scratch/user/uqpha1/logs/%x_%j_%a.error

# =============================================================================
# BAM 2X DOWNSAMPLING - ARRAY JOB VERSION
# =============================================================================
# Usage:
#   sbatch --array=0-<N-1> bam_to_2x_downsampling.sh <accession_list.txt> <bam_workdir>
#
# The accession list should contain one accession per line. Each array task finds
# <accession>.bam in <bam_workdir>, estimates depth over Chr01-Chr17, and writes
# <accession>_2x.bam back to the same directory.
# =============================================================================

set -euo pipefail

start_ts=$(date +%s)

TARGET_DEPTH=2
SEED=67
CHROM_REGEX='^Chr(0[1-9]|1[0-7])$'
THREADS="${SLURM_CPUS_PER_TASK:-4}"

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <accession_list.txt> <bam_workdir>"
    exit 1
fi

accession_list="$1"
bam_workdir="${2%/}"

if [ ! -f "$accession_list" ]; then
    echo "ERROR: Accession list not found: $accession_list"
    exit 1
fi

if [ ! -d "$bam_workdir" ]; then
    echo "ERROR: BAM workdir not found: $bam_workdir"
    exit 1
fi

if [ -z "${SLURM_ARRAY_TASK_ID:-}" ]; then
    echo "ERROR: SLURM_ARRAY_TASK_ID is not set. Submit with: sbatch --array=0-<N-1> $0 <accession_list.txt> <bam_workdir>"
    exit 1
fi

echo "SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID"
SAMPLE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$accession_list")
SAMPLE="${SAMPLE//$'\r'/}"

if [ -z "$SAMPLE" ]; then
    echo "ERROR: No accession found for task ID ${SLURM_ARRAY_TASK_ID} (line $((SLURM_ARRAY_TASK_ID + 1)))"
    exit 1
fi

module purge
module load samtools/1.18-gcc-12.3.0

bam_file="${bam_workdir}/${SAMPLE}.bam"
output_bam="${bam_workdir}/${SAMPLE}_2x.bam"

echo "=========================================="
echo "Processing accession: $SAMPLE"
echo "Input BAM: $bam_file"
echo "Output BAM: $output_bam"
echo "Target depth: ${TARGET_DEPTH}x"
echo "Threads: $THREADS"
echo "=========================================="

if [ ! -f "$bam_file" ]; then
    echo "ERROR: BAM file not found: $bam_file"
    exit 1
fi

if [ ! -f "${bam_file}.bai" ]; then
    echo "BAM index not found. Creating index: ${bam_file}.bai"
    samtools index -@ "$THREADS" "$bam_file"
fi

echo "Estimating length-weighted mean depth over Chr01-Chr17..."
depth=$(
    samtools coverage "$bam_file" | awk -v chr_re="$CHROM_REGEX" '
        $1 ~ chr_re {
            len = $3 - $2 + 1
            depth_sum += $7 * len
            length_sum += len
        }
        END {
            if (length_sum > 0) {
                printf "%.6f", depth_sum / length_sum
            } else {
                print "NA"
            }
        }'
)

if [[ ! "$depth" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
    echo "ERROR: Could not calculate depth over Chr01-Chr17 for $bam_file (depth=$depth)"
    exit 1
fi

echo "Observed depth across Chr01-Chr17: ${depth}x"

if awk -v depth="$depth" -v target="$TARGET_DEPTH" 'BEGIN { exit !(depth <= target) }'; then
    echo "ERROR: $SAMPLE is already at or below ${TARGET_DEPTH}x (observed ${depth}x). Not downsampling."
    exit 1
fi

fraction_arg=$(awk -v seed="$SEED" -v target="$TARGET_DEPTH" -v depth="$depth" 'BEGIN { printf "%.6f", seed + target / depth }')
keep_fraction=$(awk -v target="$TARGET_DEPTH" -v depth="$depth" 'BEGIN { printf "%.6f", target / depth }')

echo "Downsampling keep fraction: $keep_fraction"
echo "samtools view -s argument: $fraction_arg"

samtools view -@ "$THREADS" -b -s "$fraction_arg" "$bam_file" -o "$output_bam"
samtools index -@ "$THREADS" "$output_bam"

echo "Downsampling complete for $SAMPLE"
echo "Output BAM: $output_bam"
echo "Output index: ${output_bam}.bai"

end_ts=$(date +%s)
echo "Runtime (s): $((end_ts - start_ts))"