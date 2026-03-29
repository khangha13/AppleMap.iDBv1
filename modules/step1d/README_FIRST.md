# VCF QC + PCA (Step 1d) -- Report-Data Mode

Step1D runs a single "compute everything, decide later" workflow. The cluster
extracts all QC metrics, computes all-pairwise KING kinship, and runs PCA on
all samples. Results are exported as a self-contained Parquet report package
for local Quarto rendering. No plots, no threshold filtering, no sample
removal on the cluster.

## What It Does

1. Extract per-chromosome site metrics (biallelic SNPs only)
2. Prepare a combined VCF for PCA (`combined_for_pca.vcf.gz`)
3. PLINK2 pipeline: import, QC, KING (all pairwise), LD prune, PCA
4. Export Parquet report package + `manifest.json`
5. Tar the package for transfer

## Requirements

- `bcftools`, `plink2`, `Rscript` with `data.table`, `arrow`, `jsonlite`
- On HPC the wrappers load `miniforge/25.3.0-3`, `bcftools`, `plink/2.00a3.6-gcc-11.3.0`

## Fast Start

```bash
cd /path/to/GATK_Pipeline_KH_v1

# Interactive (prompts for VCF dir)
bash wrappers/interactive/step1d_interactive.sh

# SLURM
bash wrappers/sbatch/step1d_submit.sh [<dataset>] <vcf_dir>

# Beagle-imputed VCFs
bash wrappers/sbatch/step1d_submit.sh [<dataset>] <vcf_dir> --beagle

# Dry-run (preview)
bash wrappers/sbatch/step1d_submit.sh [<dataset>] <vcf_dir> --dry-run
```

### Deprecated Flags

The following flags are accepted with a deprecation warning but are ignored.
Step1D now runs a single combined workflow by default:

- `--qc`, `--PCA`, `--pca`, `--duplicate-check`, `--remove-relatives`

### Prepare Combined VCF (utility)

If your VCFs need cleaning (removing Chr00, standardizing contig names to
ChrNN), the prep utility runs automatically. You can also invoke it manually:

```bash
bash modules/step1d/bin/prepare_combined_for_pca.sh /path/to/vcf/directory --dry-run
bash modules/step1d/bin/prepare_combined_for_pca.sh /path/to/vcf/directory
```

## Run Options

**Interactive wrapper** (`wrappers/interactive/step1d_interactive.sh`)
Flags: `--dir=/path/to/vcfs`, `--vcf=Chr00,Chr01`, `--beagle`, `--dry-run`.

**SLURM wrapper** (`wrappers/sbatch/step1d_submit.sh`)
`bash wrappers/sbatch/step1d_submit.sh [<dataset>] <vcf_dir> [--beagle] [--dry-run]`

**Direct script call** (maximum control):
```bash
export VCF_DIR="/absolute/path/to/vcfs"
export WORK_DIR="${VCF_DIR}"
bash modules/step1d/templates/master_vcf_analysis.sh [--beagle] [--dry-run]
```

## Inputs and Patterns

- Default: 18 files `Chr00.vcf.gz` ... `Chr17.vcf.gz`. Beagle mode expects 17.
- Custom names: pass `--vcf=Chr00,Chr01` to the wrapper or set
  `VCF_INCLUDE_FILENAMES="file1.vcf.gz file2.vcf.gz"`.
- Custom pattern: set `VCF_PATTERN="chr%02d.filtered.vcf.gz"`.

## Configuration Variables

### Cache and Package Directories

```bash
STEP1D_CACHE_DIR   # Intermediate outputs (default: ${WORK_DIR}/${DATASET}_step1d_cache)
STEP1D_PACKAGE_DIR # Final Parquet package (default: ${WORK_DIR}/${DATASET}_step1d_report_package)
STEP1D_PARQUET_COMPRESSION  # Parquet compression (default: snappy)
```

### PCA Settings

```bash
STEP1D_PCA_DIR                  # PCA subdir name (default: pca_analysis)
STEP1D_PCA_MERGED_PATTERN       # Merged VCF detection (default: *merged*.vcf.gz,*merge*.vcf.gz)
STEP1D_PCA_FORCE_CONCAT         # Always concat per-chrom VCFs (default: false)
STEP1D_PCA_MERGED_EXCLUDE_CHR   # Exclude names with 'Chr' from merge detection (default: true)
```

## Output Layout

### Cache Directory (intermediate, HPC scratch)

```
${DATASET}_step1d_cache/
  Chr01_snps.vcf.gz ... Chr17_snps.vcf.gz     # Filtered SNP VCFs
  site_metrics_per_chromosome/                 # Per-chrom metric TSVs
  variant_site_metrics.tsv                     # Combined metrics
  SNP_site_meanDP.tsv                          # Mean depth (non-Beagle)
  combined_for_pca.vcf.gz                      # Combined VCF
  combined_for_pca.stats.txt                   # bcftools stats
  pca_analysis/                                # PLINK2 intermediates
    all_chromosomes.{pgen,pvar,psam}
    qc.{pgen,pvar,psam}
    king_pairwise.kin0                         # ALL pairwise KING values
    qc_pruned.{pgen,pvar,psam}
    pca.eigenvec
    pca.eigenval
```

### Report Package (final output, transfer to local)

```
${DATASET}_step1d_report_package/
  manifest.json
  qc/
    site_metrics/
      chrom=Chr01/part-000.parquet
      chrom=Chr02/part-000.parquet
      ... chrom=Chr17/
  pca/
    scores.parquet              # All samples, PC1-PC10
    variance.parquet            # 10 rows, eigenvalue + % explained
    sample_annotations.parquet  # All imported samples, passed_qc flag
    king_pairwise.parquet       # ALL N*(N-1)/2 pairs, no threshold
```

The package is also available as `${DATASET}_step1d_report_package.tar.gz`.

## Time and Resources (rough)

| Step | Small (<50 samples) | Medium (50-100) | Large (>100) |
|------|---------------------|-----------------|--------------|
| TSV extraction | 1-4 h | 4-10 h | 12-24 h |
| Combined VCF prep | minutes | minutes | 10-30 m |
| PLINK2 PCA pipeline | 1-3 h | 3-6 h | 6-12+ h |
| Parquet export | seconds | seconds | minutes |

## Troubleshooting

| Problem | Fix |
|---------|-----|
| "No VCF files specified" | Check `VCF_DIR`, naming pattern, or `VCF_INCLUDE_FILENAMES` |
| `plink2` or `bcftools` not found | `module load plink/2.00a3.6` / `module load bcftools` |
| Missing R packages | `conda activate rplot && Rscript -e "install.packages(c('data.table','arrow','jsonlite'))"` |

## File Map

- `modules/step1d/templates/master_vcf_analysis.sh` -- main orchestrator
- `modules/step1d/templates/plink2_PCA.sh` -- PCA helper (5 params: vcf_dir, rscripts_dir, plink2_bin, bcftools_bin, cache_pca_dir)
- `modules/step1d/Rscripts/export_parquet_package.R` -- Parquet exporter
- `modules/step1d/bin/prepare_combined_for_pca.sh` -- VCF preparation utility
- `wrappers/interactive/step1d_interactive.sh` -- interactive launcher
- `wrappers/sbatch/step1d_submit.sh` -- SLURM submission wrapper

**Pro tip:** dry-run before long jobs: `bash master_vcf_analysis.sh --dry-run`
