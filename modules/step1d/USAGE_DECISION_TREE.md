# Step1D Usage Decision Tree (v2 -- report-data)

**Quick Guide:** Which script should I use?

---

## Start Here

Step1D now runs a single "report-data" workflow by default. No mode flags
needed -- all computations (QC metrics, KING, PCA) run automatically.

```
Where do you want to run?
          |
    +-----+------+
    |             |
 On HPC       Locally
 (SLURM)
    |             |
run_step1d.sh   step1d_interactive.sh
```

---

## Quick Command Reference

### Scenario 1: Run on HPC (SLURM)

```bash
bash modules/step1d/bin/run_step1d.sh /path/to/vcf/directory

# Or with dataset name:
bash modules/step1d/bin/run_step1d.sh DATASET_NAME /path/to/vcf/directory
```

**What you get:**
- Site metrics (per-chrom + combined TSV)
- Combined VCF for PCA
- PLINK2 files (import, QC, KING, LD-pruned, PCA)
- Parquet report package tarball

---

### Scenario 2: Interactive (Local)

```bash
bash wrappers/interactive/step1d_interactive.sh

# Skip prompts:
bash wrappers/interactive/step1d_interactive.sh --dir=/path/to/vcfs
```

---

### Scenario 3: Beagle-Imputed Data

```bash
bash modules/step1d/bin/run_step1d.sh /path/to/vcfs --beagle

# Or interactive:
bash wrappers/interactive/step1d_interactive.sh --dir=/path/to/vcfs --beagle
```

---

### Scenario 4: Just Prepare Combined VCF

```bash
bash modules/step1d/bin/prepare_combined_for_pca.sh /path/to/vcf/directory

# Preview:
bash modules/step1d/bin/prepare_combined_for_pca.sh /path/to/vcf/directory --dry-run
```

**What you get:**
- `combined_for_pca.vcf.gz`
- `combined_for_pca.stats.txt`

---

### Scenario 5: Advanced -- Custom Environment

```bash
export VCF_DIR=/custom/path
export WORK_DIR=/output/path
export STEP1D_CACHE_DIR=/scratch/cache
export STEP1D_PACKAGE_DIR=/scratch/package

bash modules/step1d/templates/master_vcf_analysis.sh [--beagle] [--dry-run]
```

---

## Valid Flag Combinations

| Flags | What it does |
|-------|-------------|
| (none) | Full report-data workflow (default) |
| `--beagle` | Use Beagle metrics (AF, DR2, IMP) |
| `--dry-run` | Preview (works with any combination) |
| `--beagle --dry-run` | Preview Beagle mode |

### Deprecated Flags (accepted with warning)

| Flag | Old Behavior | Now |
|------|-------------|-----|
| `--qc` | QC metrics + plots only | Ignored (report-data runs all) |
| `--PCA`, `--pca` | PCA only | Ignored (report-data runs all) |
| `--duplicate-check` | KING only | Ignored (report-data runs all) |
| `--remove-relatives` | Remove relatives before PCA | Ignored (all samples in PCA) |

---

## Common Mistakes

### Mistake 1: Calling master_vcf_analysis.sh with directory argument

Wrong:
```bash
bash templates/master_vcf_analysis.sh /path/to/vcfs
```

Right:
```bash
export VCF_DIR=/path/to/vcfs
bash templates/master_vcf_analysis.sh
```

Or use a wrapper:
```bash
bash bin/run_step1d.sh /path/to/vcfs
```

---

### Mistake 2: Trying old mode flags

Wrong:
```bash
run_step1d.sh /vcfs --qc --PCA    # Error: these were mutually exclusive anyway
```

Right:
```bash
run_step1d.sh /vcfs               # report-data does everything
```

---

## Output Files

### Report Package (transfer to local machine)

```
${DATASET}_step1d_report_package/
  manifest.json
  qc/
    site_metrics/
      chrom=Chr01/part-000.parquet
      ...
  pca/
    scores.parquet
    variance.parquet
    sample_annotations.parquet
    king_pairwise.parquet
```

Also available as `${DATASET}_step1d_report_package.tar.gz`.

### Cache Directory (HPC scratch, intermediate)

```
${DATASET}_step1d_cache/
  variant_site_metrics.tsv
  combined_for_pca.vcf.gz
  pca_analysis/
    pca.eigenvec, pca.eigenval
    king_pairwise.kin0
    ...
```

---

## Environment Variables Quick Reference

### Must Set (If Not Using Wrappers)

```bash
VCF_DIR=/path/to/vcfs              # Input directory (read-only)
```

### Commonly Customized

```bash
WORK_DIR=/path/to/output           # Where to write results (default: VCF_DIR)
VCF_PATTERN="Chr%02d.vcf.gz"       # Filename pattern
STEP1D_CACHE_DIR=/scratch/cache    # Cache location
STEP1D_PACKAGE_DIR=/scratch/pkg    # Package location
```

### Advanced

```bash
STEP1D_PCA_MERGED_PATTERN="*merged*.vcf.gz"    # Merged VCF detection
STEP1D_PARQUET_COMPRESSION="snappy"             # Parquet compression
```

---

## Cheat Sheet

```bash
# Most common commands:

# 1. Run full pipeline on cluster
run_step1d.sh /vcfs

# 2. Interactive
step1d_interactive.sh

# 3. Preview what would happen
run_step1d.sh /vcfs --dry-run

# 4. Beagle-imputed data
run_step1d.sh /vcfs --beagle

# 5. Just prepare combined VCF
prepare_combined_for_pca.sh /vcfs
```

---

**Pro tip:** Start with `--dry-run` to see what would happen before running expensive analyses!
