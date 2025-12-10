#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10G
#SBATCH --job-name=vcfstat
#SBATCH --time=2-00:00:00
#SBATCH --partition=general
#SBATCH --account=a_qaafi_chs
#SBATCH -o slurm.statvcf
#SBATCH -e slurm.statvcf

# Load required modules (adjust to your cluster environment)
module purge
module load bcftools/1.18-gcc-12.3.0
module load python/3.12.3-gcccore-13.3.0
module load matplotlib/3.5.2-foss-2022a

# Input arguments
#   $1 = VCF path (required)
#   $2 = output prefix (optional; defaults to basename without .vcf/.vcf.gz)
if [[ $# -lt 1 ]]; then
  echo "Usage: bash bcfstats.sh <vcf> [prefix]" >&2
  exit 1
fi
VCF="$1"
if [[ $# -ge 2 ]]; then
  PREFIX="$2"
else
  base="$(basename "${VCF}")"
  PREFIX="${base%.vcf.gz}"
  PREFIX="${PREFIX%.vcf}"
fi

echo "Starting QC job on $VCF"

# Ensure index exists for faster random access
if [ ! -f "${VCF}.tbi" ]; then
  bcftools index -t -f "$VCF"
fi

# Count SNPs using bcftools index count
bcftools view -v snps "$VCF" -Ou | bcftools index -n > "${PREFIX}_snp_count.txt"

# Basic stats (includes depth & missingness summaries)
bcftools stats "$VCF" > "${PREFIX}_bcftools_stats.txt"

# Extract depth matrix: CHROM POS INFO/DP and per-sample DP
bcftools query -f '%CHROM\t%POS\t%INFO/DP[\t%DP]\n' "$VCF" > "${PREFIX}_depth.tsv"

# Extract per-sample missingness from GT field
bcftools query -f '[%SAMPLE\t%GT\n]' "$VCF" \
  | awk -F '\t' '{
      sample=$1; gt=$2;
      total[sample]++; if(gt==".") miss[sample]++;
    }
    END{for(s in total) print s, (miss[s]+0)/total[s];}' \
  > "${PREFIX}_missingness_individual.txt"

# Plotting with Python
python3 <<EOF
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style="whitegrid")

prefix = "${PREFIX}"

# Depth per sample (mean across sites)
depth = pd.read_csv(f"{prefix}_depth.tsv", sep="\t", header=None)
# columns: CHROM, POS, INFO_DP, sample1_DP, sample2_DP, ...
sample_cols = depth.columns[3:]
mean_depth_per_sample = depth[sample_cols].mean()

plt.figure(figsize=(8,6))
sns.histplot(mean_depth_per_sample, bins=50, color="steelblue")
plt.xlabel("Mean Depth per Individual"); plt.ylabel("Count")
plt.title("Distribution of Mean Depth per Individual")
plt.tight_layout(); plt.savefig(f"{prefix}_depth_individual_hist.png")

# Depth per site (mean across samples)
mean_depth_per_site = depth[sample_cols].mean(axis=1)
plt.figure(figsize=(8,6))
sns.histplot(mean_depth_per_site, bins=50, color="darkgreen")
plt.xlabel("Mean Depth per Site"); plt.ylabel("Count")
plt.title("Distribution of Mean Depth per Site")
plt.tight_layout(); plt.savefig(f"{prefix}_depth_site_hist.png")

# Missingness per individual (from awk output)
imiss = pd.read_csv(f"{prefix}_missingness_individual.txt", sep=" ", header=None, names=["sample","f_miss"])
plt.figure(figsize=(8,6))
sns.histplot(imiss["f_miss"], bins=50, color="firebrick")
plt.xlabel("Fraction Missing per Individual"); plt.ylabel("Count")
plt.title("Distribution of Missingness per Individual")
plt.tight_layout(); plt.savefig(f"{prefix}_missingness_individual_hist.png")
EOF
