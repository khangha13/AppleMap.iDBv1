# VCF QC + PCA (Step 1d) — Single README

This file now consolidates the old Quick Reference, How to Run PCA, and Direct Call guides into one place. For cluster onboarding and storage tips, keep `HPC_SETUP_GUIDE.md` handy.

## What It Does
- QC pipeline: extract mean depth, make per-chromosome depth and missingness plots, and combine results (18 chromosomes by default).
- PCA/duplicate workflows: merge or concat VCFs, run light QC and LD pruning, optional KING duplicate detection and relative removal, then plink2 PCA + R plots.

## Requirements
- `bcftools`, `plink2`, `Rscript` with `ggplot2`, `data.table`, `ragg`, `scales`, `ggrepel` (optional).
- On HPC the wrappers load `miniforge/25.3.0-3`, `bcftools`, `plink/2.00a3.6-gcc-11.3.0`; adjust module names as needed.

## Fast Start
### QC (default)
```bash
cd /path/to/GATK_Pipeline_KH_v1
bash wrappers/interactive/step1d_interactive.sh               # prompts for VCF dir
# or SLURM array (recommended)
bash wrappers/sbatch/step1d_submit.sh [<dataset>] <vcf_dir>
```

### PCA-only or duplicate-check
```bash
cd /path/to/GATK_Pipeline_KH_v1
bash wrappers/interactive/step1d_interactive.sh --PCA [--remove-relatives]
bash wrappers/interactive/step1d_interactive.sh --duplicate-check          # KING only

# SLURM
bash wrappers/sbatch/step1d_submit.sh [<dataset>] <vcf_dir> --PCA [--remove-relatives]
bash wrappers/sbatch/step1d_submit.sh [<dataset>] <vcf_dir> --duplicate-check
```

### Prepare Combined VCF (utility)
If your VCFs need cleaning (removing Chr00, standardizing contig names to ChrNN), use the prep utility **before** running PCA:
```bash
cd /path/to/GATK_Pipeline_KH_v1
# Dry-run to preview what would be done
bash modules/step1d/bin/prepare_combined_for_pca.sh /path/to/vcf/directory --dry-run

# Create cleaned combined_for_pca.vcf.gz (removes Chr00, standardizes to Chr01..Chr17)
bash modules/step1d/bin/prepare_combined_for_pca.sh /path/to/vcf/directory

# Then run PCA on the cleaned output
bash modules/step1d/templates/master_vcf_analysis.sh --PCA
```
**What it does:** Detects merged VCF or concatenates per-chromosome VCFs, removes Chr00 records, standardizes contigs to ChrNN format, sorts and indexes. See `--help` for details.

## Run Modes (details)
**Interactive wrapper** (`wrappers/interactive/step1d_interactive.sh`)  
Flags: `--dir=/path/to/vcfs`, `--vcf=Chr00,Chr01`, `--beagle`, `--qc` (default), `--PCA`, `--duplicate-check`, `--remove-relatives`, `--dry-run`.

**SLURM wrapper** (`wrappers/sbatch/step1d_submit.sh`)  
`bash wrappers/sbatch/step1d_submit.sh [<dataset>] <vcf_dir> [--PCA|--duplicate-check|--qc] [--beagle] [--remove-relatives]`

**Direct script call** (maximum control)  
```bash
cd /path/to/GATK_Pipeline_KH_v1
export VCF_DIR="/absolute/path/to/vcfs"
export WORK_DIR="${VCF_DIR}"
# choose one:
# export VCF_INCLUDE_FILENAMES="Chr00.vcf.gz Chr01.vcf.gz"   # explicit list
# export VCF_PATTERN="chr%02d.vcf.gz"                        # formatted pattern
# optional PCA tweaks:
export STEP1D_REMOVE_RELATIVES=true      # KING cutoff 0.125
export STEP1D_PCA_FORCE_CONCAT=true      # always concat per-chrom VCFs
bash modules/step1d/templates/master_vcf_analysis.sh --PCA   # or --qc / --duplicate-check
```
One-liner example:  
`export VCF_DIR="/data/vcfs" && export VCF_INCLUDE_FILENAMES="my.vcf.gz" && export WORK_DIR="${VCF_DIR}" && bash modules/step1d/templates/master_vcf_analysis.sh --PCA`

## Inputs & Patterns
- Default expectation: 18 files `Chr00.vcf.gz` … `Chr17.vcf.gz` (`VCF_PATTERN="Chr%02d.vcf.gz"`). Beagle mode expects 17 files (`Chr01`–`Chr17`).
- Custom names: pass `--vcf=Chr00,Chr01` to the wrapper or set `VCF_INCLUDE_FILENAMES="file1.vcf.gz file2.vcf.gz"`.
- Custom pattern: set `VCF_PATTERN="chr%02d.filtered.vcf.gz"`.
- Merged vs per-chrom VCFs (PCA): the script auto-detects merged files matching `*merged*.vcf.gz,*merge*.vcf.gz` (ignores names containing “Chr” unless `STEP1D_PCA_MERGED_EXCLUDE_CHR=false`). Set `STEP1D_PCA_FORCE_CONCAT=true` to always `bcftools concat` per-chrom VCFs.
- Output location defaults to `WORK_DIR`; PCA results land in `${WORK_DIR}/${STEP1D_PCA_DIR:-pca_analysis}`.

## PCA & Kinship Toggles
- `STEP1D_DUPLICATE_MODE=flag|remove|off` (default `flag`), `STEP1D_DUPLICATE_KING_THRESHOLD=0.45`
- `STEP1D_REMOVE_RELATIVES=true` to drop KING >0.125 before PCA (requires `--PCA`)
- `STEP1D_PCA_SHOW_LABELS=true|false`, `STEP1D_PCA_LABEL_SIZE=1.5`, `STEP1D_PCA_USE_GGREPEL=true|false`
- `STEP1D_PCA_FORCE_CONCAT=true` and `STEP1D_PCA_MERGED_EXCLUDE_CHR=true|false` control merge/concat behavior
- `STEP1D_PCA_DIR` sets the PCA output folder name

## Outputs
**QC:**  
`SNP_site_meanDP.tsv`, `depth_pdfs/Chr*.pdf`, `Chr*_missingness_vs_depth.png`, `all_chromosomes_missingness_vs_depth.png` (plus heterozygosity/QD/call rate plots when enabled).

**PCA / duplicate-check:**  
`pca_analysis/` (or `${STEP1D_PCA_DIR}/`): `pca.eigenvec`, `pca.eigenval`, `pca.eigenvec.var`, `pca_PC1_PC2.png`, `pca_PC1_PC3.png`, `pca_PC2_PC3.png`, `pca_scree.png`, and optional `king_duplicate_pairs.tsv`, `king_duplicate_samples.tsv`.

## Time & Resources (rough)
| Step | Small (<50 samples) | Medium (50–100) | Large (>100) |
|------|---------------------|-----------------|--------------|
| TSV extraction | 1–4 h | 4–10 h | 12–24 h |
| Depth plots | 5–30 m | 30–90 m | 1–4 h |
| Missingness per chr | 1–3 h | 3–6 h | 6–12 h |
| Combine plots | minutes | minutes | minutes |
| PCA | 1–3 h | 3–6 h | 6–12+ h |

## Troubleshooting (quick fixes)
| Problem | Fix |
|---------|-----|
| “No VCF files specified / VCF not found” | Check `VCF_DIR`, naming pattern, or `VCF_INCLUDE_FILENAMES`; `ls -lh "${VCF_DIR}"` |
| `plink2` or `bcftools` not found | `module load plink/2.00a3.6-gcc-11.3.0` / `module load bcftools` or set `PLINK2_BIN`/`BCFTOOLS_BIN` |
| Missing R packages | `conda activate rplot` (or your env) then install: `Rscript -e "install.packages(c('ggplot2','data.table','ragg','scales','ggrepel'), repos='https://cloud.r-project.org')"` |
| Permissions | `chmod +x *.sh` |
| Need to rerun specific chromosomes | Use array ranges: `sbatch --array=3,7,15 job_missingness_array.sh` |

## Handy CLI Snippets
```bash
# Edit / inspect config
nano vcf_analysis_config.sh
source vcf_analysis_config.sh && show_config
source vcf_analysis_config.sh && validate_config

# Monitor jobs
squeue -u $USER
seff JOB_ID
```

## File Map (most-used)
- `modules/step1d/templates/master_vcf_analysis_array.sh` — SLURM array driver (QC)
- `modules/step1d/templates/master_vcf_analysis.sh` — single-job driver (QC/PCA/duplicate)
- `modules/step1d/templates/plink2_PCA.sh` — PCA + kinship helper
- `wrappers/interactive/step1d_interactive.sh` — prompts and launches locally/HPC
- `wrappers/sbatch/step1d_submit.sh` — SLURM submission wrapper
- `modules/step1d/Rscripts/*.R` — plotting helpers

**Pro tip:** dry-run before long jobs: `bash modules/step1d/templates/master_vcf_analysis.sh --dry-run --PCA` (or `--qc`). Test with `CHR_START=0 CHR_END=1` to validate settings quickly.
