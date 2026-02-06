# Step1D Usage Decision Tree
**Quick Guide:** Which script should I use?

---

## Start Here: What do you want to do?

```
┌─────────────────────────────────────────────────────────┐
│         What do you want to analyze?                    │
└─────────────────────────────────────────────────────────┘
                           │
          ┌────────────────┴────────────────┐
          │                                 │
    ┌─────▼─────┐                     ┌────▼────┐
    │  QC plots │                     │   PCA   │
    │ + metrics │                     │analysis │
    └─────┬─────┘                     └────┬────┘
          │                                 │
          └────────────────┬────────────────┘
                           │
              ┌────────────▼───────────┐
              │  Where to run?         │
              └────────────┬───────────┘
                           │
          ┌────────────────┴────────────────┐
          │                                 │
    ┌─────▼─────┐                     ┌────▼────┐
    │  On HPC   │                     │ Locally │
    │  (SLURM)  │                     │         │
    └─────┬─────┘                     └────┬────┘
          │                                 │
          │                                 │
    Use: run_step1d.sh              Use: step1d_interactive.sh
```

---

## Quick Command Reference

### Scenario 1: Run Full QC on HPC

```bash
# Submit to SLURM, run all QC metrics + plots
bash modules/step1d/bin/run_step1d.sh /path/to/vcf/directory

# Or with dataset name:
bash modules/step1d/bin/run_step1d.sh DATASET_NAME /path/to/vcf/directory
```

**What you get:**
- Site metrics TSV
- 9 types of plots
- Filtered SNP VCFs
- Logs in SLURM output

---

### Scenario 2: Run PCA on HPC

```bash
# Submit PCA job to cluster
bash modules/step1d/bin/run_step1d.sh /path/to/vcf/directory --PCA

# With relative removal:
bash modules/step1d/bin/run_step1d.sh /path/to/vcf/directory --PCA --remove-relatives
```

**What you get:**
- `combined_for_pca.vcf.gz` (auto-prepared)
- `combined_for_pca.stats.txt` (SNP counts!)
- PLINK2 files (.pgen/.pvar/.psam)
- PCA plots (PC1 vs PC2, etc.)
- Eigenvalues/eigenvectors
- Duplicate detection results (if enabled)

---

### Scenario 3: Interactive QC (Local)

```bash
# Guided prompts for directory selection
bash wrappers/interactive/step1d_interactive.sh

# Skip prompts:
bash wrappers/interactive/step1d_interactive.sh --dir=/path/to/vcfs
```

**What you get:**
- Same as Scenario 1, but runs locally
- No SLURM submission
- Good for testing or small datasets

---

### Scenario 4: Just Prepare VCF for PCA

```bash
# Create combined_for_pca.vcf.gz without running PCA
bash modules/step1d/bin/prepare_combined_for_pca.sh /path/to/vcf/directory

# See what would happen:
bash modules/step1d/bin/prepare_combined_for_pca.sh /path/to/vcf/directory --dry-run
```

**What you get:**
- `combined_for_pca.vcf.gz`
- `combined_for_pca.stats.txt`
- Ready for downstream analysis

---

### Scenario 5: Duplicate Detection Only

```bash
# Run KING duplicate detection without PCA plots
bash modules/step1d/bin/run_step1d.sh /path/to/vcfs --duplicate-check
```

**What you get:**
- `king_duplicate_pairs.tsv` - Pairwise kinship
- `king_duplicate_samples.tsv` - Flagged sample IDs

---

### Scenario 6: Advanced - Custom Environment

```bash
# Set custom parameters
export VCF_DIR=/custom/path
export WORK_DIR=/output/path
export VCF_PATTERN="chromosome_%02d.vcf.gz"
export PLOT_IMAGE_FORMAT="pdf"
export STEP1D_PCA_MERGED_PATTERN="*combined*.vcf.gz"

# Run directly
bash modules/step1d/templates/master_vcf_analysis.sh --PCA
```

**Use case:** Non-standard VCF naming or custom output locations

---

## Flag Combinations

### Valid Combinations

| Flags | What it does |
|-------|-------------|
| `--qc` | QC metrics + plots (default) |
| `--PCA` | Combined VCF prep + PCA + plots |
| `--PCA --remove-relatives` | PCA with KING relative removal |
| `--duplicate-check` | Only duplicate detection, no PCA |
| `--beagle` | Use Beagle metrics (AF, DR2, IMP) |
| `--beagle --qc` | QC plots for imputed data |
| `--beagle --PCA` | PCA on imputed VCFs |
| `--dry-run` | Preview (works with any mode) |

### Invalid Combinations

| Flags | Error |
|-------|-------|
| `--qc --PCA` | ❌ "Multiple modes; choose one" |
| `--PCA --duplicate-check` | ❌ "Multiple modes; choose one" |
| `--qc --remove-relatives` | ❌ "--remove-relatives requires --PCA" |

---

## Common Mistakes

### Mistake 1: Calling master_vcf_analysis.sh with directory

❌ **Wrong:**
```bash
bash templates/master_vcf_analysis.sh /path/to/vcfs --PCA
```

✅ **Right:**
```bash
export VCF_DIR=/path/to/vcfs
bash templates/master_vcf_analysis.sh --PCA
```

Or use a wrapper:
```bash
bash bin/run_step1d.sh /path/to/vcfs --PCA
```

---

### Mistake 2: Mixing mode flags

❌ **Wrong:**
```bash
run_step1d.sh /vcfs --qc --PCA
```

✅ **Right:**
```bash
# Choose ONE mode:
run_step1d.sh /vcfs --qc
# OR
run_step1d.sh /vcfs --PCA
```

---

### Mistake 3: Running plink2_PCA.sh directly

❌ **Wrong:**
```bash
bash templates/plink2_PCA.sh /vcfs /rscripts
```

✅ **Right:**
```bash
# Use the wrapper:
bash bin/run_step1d.sh /vcfs --PCA
```

**Why:** plink2_PCA.sh requires 12 positional parameters in exact order!

---

### Mistake 4: Forgetting --force when re-running

❌ **Wrong:**
```bash
# Re-run after first attempt:
prepare_combined_for_pca.sh /vcfs
# Error: "Output file already exists"
```

✅ **Right:**
```bash
prepare_combined_for_pca.sh /vcfs --force
# OR remove the file manually:
rm /vcfs/combined_for_pca.vcf.gz*
prepare_combined_for_pca.sh /vcfs
```

---

## Environment Variables Quick Ref

### Must Set (If Not Using Wrappers)

```bash
VCF_DIR=/path/to/vcfs              # Input directory
```

### Commonly Customized

```bash
WORK_DIR=/path/to/output           # Where to write results (default: VCF_DIR)
VCF_PATTERN="Chr%02d.vcf.gz"       # Filename pattern
PLOT_IMAGE_FORMAT="png"            # png|pdf|svg
```

### Advanced (PCA)

```bash
STEP1D_PCA_MERGED_PATTERN="*merged*.vcf.gz"    # Merged VCF detection
STEP1D_DUPLICATE_MODE="flag"                    # off|flag|remove
STEP1D_DUPLICATE_KING_THRESHOLD="0.485"         # Kinship threshold
STEP1D_REMOVE_RELATIVES="false"                 # Remove close relatives
```

**See:** `FUNCTION_MAP.md` §6.1 for complete list (28+ variables!)

---

## Output Files

### QC Mode

```
/path/to/vcfs/
├── variant_site_metrics.tsv               # Combined metrics
├── SNP_site_meanDP.tsv                    # Depth data
├── Chr00_snps.vcf.gz ... Chr17_snps.vcf.gz  # Filtered SNPs
├── site_metrics_per_chromosome/           # Per-chr metrics
│   ├── Chr00_metrics.tsv
│   └── ...
├── depth_plots/                           # Depth vs position
│   ├── Chr00_depth_vs_position.png
│   └── ...
├── missingness_plots/                     # Missingness vs position
├── depth_vs_missingness/                  # Scatter plots
├── site_quality_plots/                    # Quality histograms
├── heterozygosity_plots/                  # Het rates
├── quality_by_depth_plots/                # QD histograms
├── call_rate_heatmaps/                    # Call rate heatmaps
└── af_distribution_plots/                 # Allele freq histograms
```

### PCA Mode

```
/path/to/vcfs/
├── combined_for_pca.vcf.gz                # Combined VCF
├── combined_for_pca.vcf.gz.csi            # Index
├── combined_for_pca.stats.txt             # SNP counts!
└── pca_analysis/                          # PCA outputs
    ├── pca.eigenvec                       # Eigenvectors
    ├── pca.eigenval                       # Eigenvalues
    ├── pca_plot_PC1_PC2.png              # Main PCA plot
    ├── pca_plot_PC3_PC4.png              # Additional PCs
    ├── king_duplicate_pairs.tsv           # Duplicate pairs (if any)
    └── king_duplicate_samples.tsv         # Flagged samples
```

---

## Cheat Sheet

```bash
# Most common commands:

# 1. Full QC on cluster (default mode)
run_step1d.sh /vcfs

# 2. PCA on cluster
run_step1d.sh /vcfs --PCA

# 3. Interactive QC
step1d_interactive.sh

# 4. Preview what would happen
run_step1d.sh /vcfs --dry-run

# 5. Just prepare combined VCF
prepare_combined_for_pca.sh /vcfs

# 6. Check for duplicates
run_step1d.sh /vcfs --duplicate-check

# 7. Beagle-imputed data QC
run_step1d.sh /vcfs --beagle
```

---

**Pro tip:** Start with `--dry-run` to see what would happen before running expensive analyses!
