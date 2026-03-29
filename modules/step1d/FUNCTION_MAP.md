# Step1D Function and Flag Map (v2 -- report-data)

**Updated:** 2026-02-03
**Purpose:** Comprehensive mapping of all entry points and flags after the
Parquet report package refactor.

---

## Summary

Step1D now runs a single "compute everything, decide later" workflow. Old mode
flags (`--qc`, `--PCA`, `--duplicate-check`, `--remove-relatives`) are accepted
with a deprecation warning and mapped to the default report-data workflow.

The cluster extracts metrics, computes KING (all pairwise), runs PCA on all
samples, and exports a Parquet package for local Quarto rendering.

---

## 1. Entry Points Overview

| Entry Point | Type | Primary Use Case | User Facing |
|-------------|------|------------------|-------------|
| **run_step1d.sh** | SLURM wrapper | Submit job to cluster | Yes |
| **step1d_interactive.sh** | Interactive wrapper | Local/interactive analysis | Yes |
| **master_vcf_analysis.sh** | Main orchestrator | Actual computation | Advanced |
| **plink2_PCA.sh** | PCA helper | PLINK2 pipeline (5 params) | Internal |
| **prepare_combined_for_pca.sh** | VCF prep utility | Prepare combined VCF | Yes |
| **export_parquet_package.R** | Parquet exporter | Write report package | Internal |

---

## 2. Flag Matrix

| Flag | run_step1d.sh | step1d_interactive.sh | master_vcf_analysis.sh |
|------|--------------|---------------------|----------------------|
| `--help`, `-h` | via usage | Yes | Yes |
| `--dry-run`, `-n` | Yes | Yes | Yes |
| `--beagle` | Yes | Yes | Yes |
| `--dir=PATH` | N/A | Yes | N/A |
| `--vcf=LIST` | N/A | Yes | N/A |
| `--qc` | Deprecated | Deprecated | Deprecated |
| `--PCA`, `--pca` | Deprecated | Deprecated | Deprecated |
| `--duplicate-check` | Deprecated | Deprecated | Deprecated |
| `--remove-relatives` | Deprecated | Deprecated | Deprecated |

**Deprecated flags** are accepted with a warning and ignored. The default
`report-data` workflow runs all computations automatically.

---

## 3. Script Details

### 3.1 run_step1d.sh (SLURM Wrapper)

**Signature:**
```bash
run_step1d.sh <dataset_name> <vcf_directory> [--beagle] [--dry-run]
run_step1d.sh <vcf_directory> [--beagle] [--dry-run]
```

**What it does:**
1. Validates dataset/VCF directory
2. Reads Step1D config from pipeline_config.sh
3. Generates SLURM script from template
4. Submits job with sbatch
5. Returns job ID

---

### 3.2 step1d_interactive.sh (Interactive Wrapper)

**Signature:**
```bash
step1d_interactive.sh [--dir=PATH] [--vcf=LIST] [--beagle] [--dry-run]
```

**Special Features:**
- VCF auto-detection (finds Chr*.vcf.gz automatically)
- VCF filtering (`--vcf=Chr01,Chr02` to analyze subset)
- Interactive prompts for directory selection
- Fallback logic for non-standard naming

---

### 3.3 master_vcf_analysis.sh (Main Orchestrator)

**Signature:**
```bash
master_vcf_analysis.sh [--beagle] [--dry-run]
```

Does NOT accept directory arguments. Uses environment variables:
```bash
VCF_DIR=/path/to/vcfs              # REQUIRED (read-only input)
WORK_DIR=/path/to/output           # Optional, defaults to VCF_DIR
```

**Linear Workflow:**
1. Extract per-chromosome site metrics (bcftools query)
2. Assemble combined `variant_site_metrics.tsv`
3. Derive mean depth TSV (non-Beagle only)
4. Prepare `combined_for_pca.vcf.gz` via prepare_combined_for_pca.sh
5. Run plink2_PCA.sh (import, QC, KING, LD prune, PCA)
6. Run export_parquet_package.R (Parquet export)
7. Create tarball

**Beagle Mode Differences:**
- Different metrics: AF, DR2, IMP (no depth, no QD, no MQ)
- Same pipeline flow; Parquet schema fills missing columns with null

---

### 3.4 plink2_PCA.sh (PCA Helper -- Internal)

**Signature (5 positional params):**
```bash
plink2_PCA.sh <vcf_dir> <rscripts_dir> [plink2_bin] [bcftools_bin] [cache_pca_dir]
```

**Pipeline steps (each with stamp-based skip check):**
1. Import VCF to PLINK (`--snps-only just-acgt --max-alleles 2`)
2. QC filter (`--geno 0.05 --mind 0.10 --maf 0.01`)
3. KING table (`--make-king-table`, preserves raw `.kin0`)
4. LD prune (`--indep-pairwise 200 50 0.2`)
5. PCA on all pruned samples (`--pca 10 biallelic-var-wts`)

No duplicate filtering, no relative removal, no plotting.

---

### 3.5 export_parquet_package.R (Parquet Exporter -- Internal)

**Arguments:**
```bash
Rscript export_parquet_package.R \
  --cache-dir /path/to/cache \
  --package-dir /path/to/package \
  --dataset-name MyDataset \
  --beagle false \
  --compression snappy
```

**Dependencies:** `data.table`, `arrow`, `jsonlite`

**Exports:**
1. QC site metrics -> Hive-partitioned Parquet by chrom
2. PCA scores -> `pca/scores.parquet`
3. PCA variance -> `pca/variance.parquet`
4. Sample annotations -> `pca/sample_annotations.parquet`
5. King pairwise (all pairs) -> `pca/king_pairwise.parquet`
6. Manifest -> `manifest.json`

---

### 3.6 prepare_combined_for_pca.sh (VCF Prep Utility)

**Signature:**
```bash
prepare_combined_for_pca.sh <vcf_directory> [--force] [--dry-run] [--help]
```

**What it does:**
1. Detect merged VCF or concatenate per-chromosome VCFs
2. Remove Chr00 records
3. Standardize contigs to ChrNN format
4. Sort and index (CSI)
5. Generate bcftools stats

---

## 4. Configuration Variables

### New in v2

```bash
STEP1D_CACHE_DIR              # Cache directory (default: set by master based on WORK_DIR)
STEP1D_PACKAGE_DIR            # Package directory (default: set by master based on WORK_DIR)
STEP1D_PARQUET_COMPRESSION    # Parquet compression (default: snappy)
```

### Kept from v1

```bash
STEP1D_PCA_DIR                # PCA subdir name (default: pca_analysis)
STEP1D_PCA_MERGED_PATTERN     # Merged VCF detection pattern
STEP1D_PCA_FORCE_CONCAT       # Always concat per-chrom VCFs
STEP1D_PCA_MERGED_EXCLUDE_CHR # Exclude names with Chr from merge detection
PLINK2_BIN                    # plink2 binary path
BCFTOOLS_BIN                  # bcftools binary path
```

### Removed in v2

The following config variables have been removed:
- `STEP1D_DUPLICATE_MODE` (was: flag/remove/off)
- `STEP1D_DUPLICATE_KING_THRESHOLD` (was: 0.45/0.485)
- `STEP1D_PCA_ONLY` (was: false)
- `STEP1D_REMOVE_RELATIVES` (was: false)
- `STEP1D_PCA_SHOW_LABELS` (was: true)
- `STEP1D_PCA_LABEL_SIZE` (was: 1.5)
- `STEP1D_PCA_USE_GGREPEL` (was: true)
- `STEP1D_AF_PLOTS_DIR` (was: af_distribution_plots)
- `STEP1D_AF_HIST_BINS` (was: 50)

All filtering and visualization decisions now happen on the local Quarto side.

---

## 5. What to Use When

| Scenario | Command |
|----------|---------|
| Report data on cluster (SLURM) | `run_step1d.sh /vcfs` |
| Interactive | `step1d_interactive.sh` |
| Beagle data | `run_step1d.sh /vcfs --beagle` |
| Preview | `run_step1d.sh /vcfs --dry-run` |
| Prepare combined VCF only | `prepare_combined_for_pca.sh /vcfs` |
| Advanced/custom | `export VCF_DIR=/vcfs; master_vcf_analysis.sh` |

---

## 6. Changes from v1

- **Removed:** All plot generation (8 R plotting scripts no longer called)
- **Removed:** `ggplot2`, `ragg`, `scales`, `ggrepel` R dependencies
- **Removed:** Duplicate filtering/removal and relative removal from cluster
- **Removed:** Mode branching (`--qc`/`--PCA`/`--duplicate-check`)
- **Added:** Parquet export via `export_parquet_package.R`
- **Added:** `data.table`, `arrow`, `jsonlite` R dependencies
- **Added:** `STEP1D_CACHE_DIR` and `STEP1D_PACKAGE_DIR`
- **Simplified:** `plink2_PCA.sh` from 12 positional params to 5
- **Simplified:** `master_vcf_analysis.sh` from multi-mode to single linear flow

---

**Report compiled by:** AI Assistant
**Version:** 2.0 (report-data refactor)
