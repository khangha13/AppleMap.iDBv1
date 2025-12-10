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
#Execute
module purge
module load vcftools/0.1.16-gcc-11.3.0
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

# Run the QC pipeline
echo "Starting QC job on $VCF"

# Count SNPs
bcftools view -v snps $VCF | grep -v "^#" | wc -l > ${PREFIX}_snp_count.txt

# Depth calculations
vcftools --gzvcf $VCF --depth --out ${PREFIX}
vcftools --gzvcf $VCF --site-mean-depth --out ${PREFIX}
vcftools --gzvcf $VCF --site-depth --out ${PREFIX}

# Missingness calculations
vcftools --gzvcf $VCF --missing-indv --out ${PREFIX}
vcftools --gzvcf $VCF --missing-site --out ${PREFIX}

# Stats summary
bcftools stats $VCF > ${PREFIX}_bcftools_stats.txt

# Plotting with Python
python3 <<EOF
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style="whitegrid")

depth = pd.read_csv("${PREFIX}.depth.mean", delim_whitespace=True)
plt.figure(figsize=(8,6))
sns.histplot(depth['MEAN_DEPTH'], bins=50, color="steelblue")
plt.xlabel("Mean Depth per Individual"); plt.ylabel("Count")
plt.title("Distribution of Mean Depth per Individual")
plt.tight_layout(); plt.savefig("${PREFIX}_depth_individual_hist.png")

site_depth = pd.read_csv("${PREFIX}.site.depth.mean", delim_whitespace=True)
plt.figure(figsize=(8,6))
sns.histplot(site_depth['MEAN_DEPTH'], bins=50, color="darkgreen")
plt.xlabel("Mean Depth per Site"); plt.ylabel("Count")
plt.title("Distribution of Mean Depth per Site")
plt.tight_layout(); plt.savefig("${PREFIX}_depth_site_hist.png")

imiss = pd.read_csv("${PREFIX}.imiss", delim_whitespace=True)
plt.figure(figsize=(8,6))
sns.histplot(imiss['F_MISS'], bins=50, color="firebrick")
plt.xlabel("Fraction Missing per Individual"); plt.ylabel("Count")
plt.title("Distribution of Missingness per Individual")
plt.tight_layout(); plt.savefig("${PREFIX}_missingness_individual_hist.png")

lmiss = pd.read_csv("${PREFIX}.lmiss", delim_whitespace=True)
plt.figure(figsize=(8,6))
sns.histplot(lmiss['F_MISS'], bins=50, color="purple")
plt.xlabel("Fraction Missing per Site"); plt.ylabel("Count")
plt.title("Distribution of Missingness per Site")
plt.tight_layout(); plt.savefig("${PREFIX}_missingness_site_hist.png")
EOF

echo "QC job finished. Results saved with prefix: $PREFIX"