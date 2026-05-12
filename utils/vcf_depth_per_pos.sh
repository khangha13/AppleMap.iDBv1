#!/bin/bash --login
# =============================================================================
# SLURM (optional): submit with working directory = where you want logs, e.g.
#   sbatch vcf_depth_per_pos.sh /path/to/vcfs [/path/to/output]
# Resources align with Step 1A SLURM defaults in config/pipeline_config.sh
# (CPUS/MEMORY/account/partition); walltime matches other VCF utils here (bcfstats).
# Override at submit time: sbatch --mem=64G --time=7-00:00:00 vcf_depth_per_pos.sh ...
# =============================================================================
#SBATCH --job-name=vcf_snp_depth
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=60G
#SBATCH --time=5-00:00:00
#SBATCH --partition=general
#SBATCH --account=a_qaafi_cas
#SBATCH -o /scratch/project/bigdata_apple/logs/snp_depth_%j.out
#SBATCH -e /scratch/project/bigdata_apple/logs/snp_depth_%j.err
# =============================================================================
# vcf_depth_per_pos.sh — Mean SNP depth + missingness per position (Chr01–Chr17)
# =============================================================================
# USAGE:
#   bash vcf_depth_per_pos.sh <vcf_dir> [output_dir] [vcf_pattern]
#   sbatch vcf_depth_per_pos.sh <vcf_dir> [output_dir] [vcf_pattern]
#
#   vcf_pattern  optional glob/substring to filter VCF filenames (default: *.vcf.gz)
#   Examples:
#     bash vcf_depth_per_pos.sh /data/vcfs                              # all .vcf.gz
#     bash vcf_depth_per_pos.sh /data/vcfs ./out *_noCommon.vcf.gz     # only noCommon
#     bash vcf_depth_per_pos.sh /data/vcfs ./out *_consolidated_noCommon.vcf.gz
#
# WHAT IT DOES:
#   1. Streams all Chr01–Chr17 VCFs once (FORMAT/DP and FORMAT/GT together)
#   2. Per SNP row: writes CHROM, POS, mean_dp (mean FORMAT/DP across samples)
#      and missingness (% samples with GT=./. at that site)
#      and accumulates per-sample sums for genome-wide mean depth
#   3. Writes snp_depth.tsv (4 cols) and snp_depth.rds (list: sites, per_sample)
#
# OUTPUT:
#   <output_dir>/snp_depth.tsv — CHROM  POS  mean_dp  missingness
#   <output_dir>/snp_depth.rds — list(sites = data.table, per_sample = data.table)
#
# NOTE: If snp_depth.tsv already exists from a previous run (3-column format),
#   delete it (and snp_depth.rds) before re-running to get the new 4-column format.
#
# REQUIREMENTS:
#   - bcftools (module or PATH)
#   - Rscript + data.table R package
# =============================================================================

set -euo pipefail
exec 2>&1   # merge stderr into stdout so all output goes to the .out log

VCF_DIR="${1:?Usage: bash $0 <vcf_dir> [output_dir] [vcf_pattern]}"
OUTPUT_DIR="${2:-${VCF_DIR}/depth}"
VCF_PATTERN="${3:-*.vcf.gz}"

mkdir -p "$OUTPUT_DIR"

TSV_OUT="$OUTPUT_DIR/snp_depth.tsv"
RDS_OUT="$OUTPUT_DIR/snp_depth.rds"
PS_TMP="$OUTPUT_DIR/.per_sample_mean_depth.tmp.tsv"

T0=$(date +%s)

# =============================================================================
# Modules
# =============================================================================
module purge
module load bcftools/1.18-gcc-12.3.0
module load miniforge/26.1.0-0
source "$ROOTMINIFORGE/etc/profile.d/conda.sh"
conda activate rplot

# =============================================================================
# Find VCFs matching the requested pattern, sorted
# =============================================================================
mapfile -t VCFS < <(ls "$VCF_DIR"/$VCF_PATTERN 2>/dev/null | sort)
if [[ ${#VCFS[@]} -eq 0 ]]; then
    echo "[$(date '+%H:%M:%S')] ERROR: no files matching '$VCF_PATTERN' found in $VCF_DIR"
    exit 1
fi

echo "[$(date '+%H:%M:%S')] Starting"
echo "[$(date '+%H:%M:%S')] Pattern : $VCF_PATTERN"
echo "[$(date '+%H:%M:%S')] Found   : ${#VCFS[@]} VCF file(s) in $VCF_DIR"
echo "[$(date '+%H:%M:%S')] Output  : $OUTPUT_DIR"

# =============================================================================
# Index any unindexed VCFs
# =============================================================================
chr_in_range() {
    local num
    num=$(basename "$1" | grep -oP '(?i)(?<=chr)\d+' | head -1)
    [[ -n "$num" ]] && (( 10#$num >= 1 && 10#$num <= 17 ))
}

for VCF in "${VCFS[@]}"; do
    chr_in_range "$VCF" || continue
    if [[ ! -f "${VCF}.csi" && ! -f "${VCF}.tbi" ]]; then
        echo "[$(date '+%H:%M:%S')] Indexing $(basename "$VCF")..."
        bcftools index -c "$VCF"
    fi
done

# First Chr01–Chr17 VCF defines sample column order (must match all chrom VCFs)
REF_VCF=""
for f in "${VCFS[@]}"; do
    if chr_in_range "$f"; then
        REF_VCF="$f"
        break
    fi
done
if [[ -z "$REF_VCF" ]]; then
    echo "[$(date '+%H:%M:%S')] ERROR: no Chr01–Chr17 .vcf.gz found in $VCF_DIR"
    exit 1
fi

SAMPLES_TAB=$(bcftools query -l "$REF_VCF" | paste -sd $'\t' -)

# Stream all SNP rows from all chromosomes — DP and GT interleaved per sample:
# $1=CHROM $2=POS  then for sample i (1-based): dp=$(2i+1)  gt=$(2i+2)
# Progress echoes go to stderr so they don't corrupt the stdout pipe to awk.
stream_snps() {
    for VCF in "${VCFS[@]}"; do
        chr_in_range "$VCF" || continue
        echo "[$(date '+%H:%M:%S')] Streaming : $(basename "$VCF")" >&2
        bcftools query -f '%CHROM\t%POS[\t%DP\t%GT]\n' -i 'TYPE="snp"' "$VCF"
        echo "[$(date '+%H:%M:%S')] Done      : $(basename "$VCF")" >&2
    done
}

# =============================================================================
# Build TSV + per-sample accumulator in one awk (single bcftools stream)
# Column layout per row: $1=CHROM $2=POS  sample i: dp=$(2i+1)  gt=$(2i+2)
# TSV output: CHROM  POS  mean_dp  missingness(%)
# =============================================================================
run_combined_awk() {
    # shellcheck disable=SC2016
    stream_snps | awk -v TSVOUT="$TSV_OUT" -v PSOUT="$PS_TMP" -v ST="$SAMPLES_TAB" '
        BEGIN {
            ns = split(ST, name, "\t")
        }
        {
            sm = 0; nn = 0; miss = 0
            for (i = 1; i <= ns; i++) {
                dp = $(2*i + 1)
                gt = $(2*i + 2)
                # missingness: GT is missing
                if (gt == "./." || gt == ".|." || gt == ".") {
                    miss++
                }
                # depth accumulation (skip missing DP)
                if (dp != "." && dp != "NA") {
                    sm += dp
                    nn++
                    sum[i] += dp
                    cnt[i]++
                }
            }
            miss_pct = miss / ns * 100
            if (nn > 0)
                printf "%s\t%s\t%.2f\t%.2f\n", $1, $2, sm / nn, miss_pct >> TSVOUT
            else
                printf "%s\t%s\tNA\t%.2f\n", $1, $2, miss_pct >> TSVOUT
        }
        END {
            for (i = 1; i <= ns; i++)
                if (cnt[i] > 0)
                    printf "%s\t%.6f\n", name[i], sum[i] / cnt[i] >> PSOUT
        }
    '
}

# Only per-sample means (when TSV already exists but RDS is missing)
run_ps_only_awk() {
    stream_snps | awk -v PSOUT="$PS_TMP" -v ST="$SAMPLES_TAB" '
        BEGIN { ns = split(ST, name, "\t") }
        {
            for (i = 1; i <= ns; i++) {
                dp = $(2*i + 1)
                if (dp != "." && dp != "NA") {
                    sum[i] += dp
                    cnt[i]++
                }
            }
        }
        END {
            for (i = 1; i <= ns; i++)
                if (cnt[i] > 0)
                    printf "%s\t%.6f\n", name[i], sum[i] / cnt[i] >> PSOUT
        }
    '
}

echo "[$(date '+%H:%M:%S')] Extracting SNP depth (one pass; TSV + per-sample means)..."

if [[ -f "$TSV_OUT" ]]; then
    echo "[$(date '+%H:%M:%S')] TSV exists — skipping site rows ($(du -sh "$TSV_OUT" | cut -f1))"
else
    printf '%b\n' "CHROM\tPOS\tmean_dp\tmissingness" > "$TSV_OUT"
    rm -f "$PS_TMP"
    T_ALL=$(date +%s)
    run_combined_awk
    echo "[$(date '+%H:%M:%S')] TSV + per-sample accum ($(($(date +%s) - T_ALL))s)"
    echo "[$(date '+%H:%M:%S')] TSV saved ($(du -sh "$TSV_OUT" | cut -f1))"
fi

# =============================================================================
# RDS — list(sites, per_sample)
# =============================================================================
if [[ -f "$RDS_OUT" ]]; then
    echo "[$(date '+%H:%M:%S')] RDS exists — skipping"
else
    # If TSV was reused from an older run, PS_TMP was never written — re-stream once.
    if [[ ! -s "$PS_TMP" ]]; then
        echo "[$(date '+%H:%M:%S')] Re-streaming VCFs for per-sample means (needed for RDS)..."
        T2=$(date +%s)
        rm -f "$PS_TMP"
        run_ps_only_awk
        echo "[$(date '+%H:%M:%S')] Per-sample pass done ($(($(date +%s) - T2))s)"
    fi

    if [[ ! -s "$PS_TMP" ]]; then
        echo "[$(date '+%H:%M:%S')] ERROR: per-sample summary empty. Check VCFs / bcftools."
        exit 1
    fi

    echo "[$(date '+%H:%M:%S')] Writing RDS..."

    Rscript -e "
      suppressPackageStartupMessages(library(data.table))
      sites <- fread('$TSV_OUT')
      per_sample <- fread('$PS_TMP', col.names = c('sample', 'mean_depth'))
      saveRDS(list(sites = sites, per_sample = per_sample), '$RDS_OUT')
    "
    rm -f "$PS_TMP"
    echo "[$(date '+%H:%M:%S')] RDS saved ($(du -sh "$RDS_OUT" | cut -f1))"
fi

rm -f "$PS_TMP"

ELAPSED=$(( $(date +%s) - T0 ))
echo "[$(date '+%H:%M:%S')] Done — ${ELAPSED}s elapsed"
echo "[$(date '+%H:%M:%S')] TSV : $TSV_OUT"
echo "[$(date '+%H:%M:%S')] RDS : $RDS_OUT"
