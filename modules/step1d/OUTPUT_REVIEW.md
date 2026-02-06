# Step1D Output Review

**Generated:** February 3, 2026  
**Purpose:** Document SLURM logging and all output files produced by Step1D

---

## 1. SLURM Output Files (.out and .err)

### Location
SLURM job logs are placed in a **dataset-specific subdirectory** under the configured log base:

```
${LOG_BASE_PATH}/${dataset_name}/slurm/
```

**Resolution Logic:**
1. If `PIPELINE_LOG_DIR_OVERRIDE` is set → use that as base
2. Else if `LOG_BASE_PATH` is set → use `${LOG_BASE_PATH}/${dataset_name}/slurm/`
3. Else default to → `${PIPELINE_SCRATCH_BASE}/logs/${dataset_name}/slurm/`

**Example with defaults:**
```bash
# For dataset "Apple_2024"
/scratch/user/uqpha1/logs/Apple_2024/slurm/
```

### File Naming Pattern
```
1D_<JOBID>_<ARRAY_ID>_<JOB_NAME>_<UNIQUE>.output
1D_<JOBID>_<ARRAY_ID>_<JOB_NAME>_<UNIQUE>.error
```

**Format:** `${step_name}_%A_%a_%x_%j.{output,error}`

Where:
- `%A` = Main job ID (e.g., `12345678`)
- `%a` = Array task ID (e.g., `0` for step1d, which uses array=0-0)
- `%x` = Job name (e.g., `Apple_GATK_1D_Apple_2024`)
- `%j` = Job allocation number

**Example:**
```
1D_12345678_0_Apple_GATK_1D_Apple_2024_12345678.output
1D_12345678_0_Apple_GATK_1D_Apple_2024_12345678.error
```

### SLURM Script Location
Generated SLURM scripts are stored in:
```
${PIPELINE_SLURM_SCRIPT_DIR}/Apple_GATK_1D_${dataset_name}_<timestamp>.sh
```

Default: `${PIPELINE_ROOT}/slurm_scripts/`

### Configuration Source
SLURM output paths are configured in:
- **lib/slurm.sh** (lines 82-104): Creates `#SBATCH` directives
- **lib/logging.sh** (lines 43-72): `resolve_log_root()` function
- **config/environment.sh**: Sets `PIPELINE_LOG_BASE` or `PIPELINE_SCRATCH_BASE`

---

## 2. Folder/Directory Names

### Overview
Step1D creates multiple output directories depending on the mode:

| Directory | Purpose | Default Name | Configurable Via |
|-----------|---------|--------------|------------------|
| Metrics by chromosome | Per-chromosome TSV metrics | `metrics_by_chrom/` | `METRICS_BY_CHROM_DIR` |
| Depth PDFs | Per-chromosome depth plots (PDF) | `depth_pdfs/` | `DEPTH_PLOT_DIR` |
| Missingness plots | Per-chromosome missingness PNG/PDF | `missingness_plots/` | `MISSINGNESS_PLOTS_DIR` |
| Depth vs Missingness | Per-chromosome depth-vs-miss plots | `depth_vs_missingness/` | `DEPTH_MISSINGNESS_DIR` |
| PCA analysis | PCA outputs and plots | `pca_analysis/` | `STEP1D_PCA_DIR` |
| SLURM logs | Job output/error files | `${LOG_BASE_PATH}/${dataset}/slurm/` | `PIPELINE_LOG_DIR_OVERRIDE` |

### Default Layout (QC Mode)
```
${WORK_DIR}/
├── metrics_by_chrom/
│   ├── Chr00_metrics.tsv
│   ├── Chr01_metrics.tsv
│   └── ...Chr17_metrics.tsv
├── depth_pdfs/
│   ├── Chr00.pdf
│   ├── Chr01.pdf
│   └── ...Chr17.pdf
├── missingness_plots/
│   ├── Chr00_missingness_vs_position.png
│   ├── Chr01_missingness_vs_position.png
│   └── ...Chr17_missingness_vs_position.png
├── depth_vs_missingness/
│   ├── Chr00_depth_vs_missingness.png
│   ├── Chr01_depth_vs_missingness.png
│   └── ...Chr17_depth_vs_missingness.png
├── SNP_site_meanDP.tsv
├── variant_site_metrics.tsv
└── all_chromosomes_missingness_vs_depth.png
```

### Default Layout (PCA Mode)
```
${WORK_DIR}/
├── combined_for_pca.vcf.gz          # If auto-prepared
├── combined_for_pca.vcf.gz.csi      # Index
├── combined_for_pca.stats.txt       # bcftools stats
└── pca_analysis/
    ├── pca.bed                      # PLINK2 binary genotypes
    ├── pca.bim                      # PLINK2 variant info
    ├── pca.fam                      # PLINK2 sample info
    ├── pca.eigenvec                 # PCA sample coordinates
    ├── pca.eigenval                 # PCA variance explained
    ├── pca.eigenvec.var             # Variance per PC
    ├── pca_PC1_PC2.png              # PC1 vs PC2 scatter plot
    ├── pca_scree.png                # Scree plot (variance explained)
    ├── king_duplicate_pairs.tsv     # Duplicate pairs (if detected)
    └── king_duplicate_samples.tsv   # Flagged duplicate samples
```

### Configuration Variables
All directory names can be customized via environment variables:

```bash
export METRICS_BY_CHROM_DIR="my_custom_metrics/"
export DEPTH_PLOT_DIR="my_depth_plots/"
export MISSINGNESS_PLOTS_DIR="my_miss_plots/"
export DEPTH_MISSINGNESS_DIR="my_depth_vs_miss/"
export STEP1D_PCA_DIR="my_pca_results/"
```

---

## 3. All Outputs Produced by Step1D

### QC Mode Outputs (`--qc`, default)

#### TSV Data Files
| File | Location | Description |
|------|----------|-------------|
| `SNP_site_meanDP.tsv` | `${WORK_DIR}/` | Mean depth per site across samples |
| `variant_site_metrics.tsv` | `${WORK_DIR}/` | Combined variant metrics (DP, QD, MQ, etc.) |
| `Chr*_metrics.tsv` | `metrics_by_chrom/` | Per-chromosome metrics (one file per chromosome) |

**Columns in variant_site_metrics.tsv:**
- CHROM, POS, ID, REF, ALT
- DP (depth), QD (quality by depth), MQ (mapping quality)
- QUAL, BaseQRankSum, MQRankSum, ReadPosRankSum
- FS (Fisher strand), SOR (strand odds ratio)
- InbreedingCoeff
- (Additional INFO fields if present)

#### Plot Files (PNG/PDF)
| File Pattern | Location | Description |
|--------------|----------|-------------|
| `Chr*.pdf` | `depth_pdfs/` | Per-chromosome depth distribution plots (PDF) |
| `Chr*_missingness_vs_position.png` | `missingness_plots/` | Per-chromosome missingness by position |
| `Chr*_depth_vs_missingness.png` | `depth_vs_missingness/` | Per-chromosome depth vs missingness scatter |
| `all_chromosomes_missingness_vs_depth.png` | `${WORK_DIR}/` | Combined plot across all chromosomes |

**Note:** Plot format can be changed via `PLOT_IMAGE_FORMAT=pdf` (default: `png`)

#### Optional QC Plots (if enabled)
- **Heterozygosity plots** (via `plot_heterozygosity.R`)
- **Quality by depth** (via `plot_quality_by_depth.R`)
- **Call rate heatmaps** (via `plot_call_rate_heatmap.R`)
- **Allele frequency distributions** (via `plot_af_distribution.R`)
- **Site quality plots** (via `plot_site_quality.R`)

---

### PCA Mode Outputs (`--PCA`)

#### VCF Preparation (auto-generated if needed)
| File | Location | Description |
|------|----------|-------------|
| `combined_for_pca.vcf.gz` | `${VCF_DIR}/` | Merged/concatenated VCF (Chr00 removed, ChrNN standardized) |
| `combined_for_pca.vcf.gz.csi` | `${VCF_DIR}/` | CSI index |
| `combined_for_pca.stats.txt` | `${VCF_DIR}/` | bcftools stats summary |

**Preparation handled by:** `modules/step1d/bin/prepare_combined_for_pca.sh`

#### PLINK2 Binary Files
| File | Location | Description |
|------|----------|-------------|
| `pca.bed` | `${STEP1D_PCA_DIR}/` | Binary genotype matrix |
| `pca.bim` | `${STEP1D_PCA_DIR}/` | Variant information (chr, ID, pos, alleles) |
| `pca.fam` | `${STEP1D_PCA_DIR}/` | Sample information (FID, IID) |
| `pca.log` | `${STEP1D_PCA_DIR}/` | PLINK2 log file |

**QC Filters Applied:**
- `--geno 0.1` (exclude variants missing in >10% samples)
- `--mind 0.1` (exclude samples missing >10% variants)
- `--maf 0.01` (exclude variants with MAF <1%)

**LD Pruning:**
- `--indep-pairwise 50 10 0.2` (50kb window, 10 variant step, r² > 0.2)

#### PCA Results
| File | Location | Description | Format |
|------|----------|-------------|--------|
| `pca.eigenvec` | `${STEP1D_PCA_DIR}/` | Sample coordinates (PC1-PC10) | TSV: FID, IID, PC1, PC2, ... |
| `pca.eigenval` | `${STEP1D_PCA_DIR}/` | Eigenvalues (variance per PC) | One value per line (PC1-PC10) |
| `pca.eigenvec.var` | `${STEP1D_PCA_DIR}/` | Variance explained per PC | Single line: 10 values |

#### PCA Plots
| File | Location | Description |
|------|----------|-------------|
| `pca_PC1_PC2.png` | `${STEP1D_PCA_DIR}/` | Scatter plot of PC1 vs PC2 (1600×1200, 200 DPI) |
| `pca_scree.png` | `${STEP1D_PCA_DIR}/` | Scree plot showing variance explained by each PC |

**Plot Features:**
- Sample labels (configurable via `STEP1D_PCA_SHOW_LABELS`)
- ggrepel for non-overlapping labels (via `STEP1D_PCA_USE_GGREPEL`)
- Duplicate highlighting (if detected)
- Customizable label size (`STEP1D_PCA_LABEL_SIZE`)

#### Duplicate Detection Outputs (if enabled)
| File | Location | Description |
|------|----------|-------------|
| `king_duplicate_pairs.tsv` | `${STEP1D_PCA_DIR}/` | All duplicate pairs (KING > threshold) |
| `king_duplicate_samples.tsv` | `${STEP1D_PCA_DIR}/` | Flagged samples to remove |
| `king_duplicates.king` | `${STEP1D_PCA_DIR}/` | Raw KING output |
| `king_duplicates.king.id` | `${STEP1D_PCA_DIR}/` | Sample IDs from KING |

**Thresholds:**
- Duplicates: KING > 0.485 (configurable via `STEP1D_DUPLICATE_KING_THRESHOLD`)
- Close relatives: KING > 0.125 (for `--remove-relatives` mode)

**Modes:**
- `flag` (default): Detect and report duplicates, include in PCA
- `remove`: Remove duplicates before PCA
- `off`: Skip duplicate detection

---

### Duplicate Check Mode Outputs (`--duplicate-check`)

Same as PCA duplicate detection outputs, but runs **only** KING analysis without PCA:

```
${STEP1D_PCA_DIR}/
├── pca.bed
├── pca.bim
├── pca.fam
├── king_duplicate_pairs.tsv
├── king_duplicate_samples.tsv
├── king_duplicates.king
└── king_duplicates.king.id
```

---

### Utility: Combined VCF Preparation Outputs

When running `prepare_combined_for_pca.sh` manually or via PCA workflow:

| File | Location | Size (typical) | Description |
|------|----------|----------------|-------------|
| `combined_for_pca.vcf.gz` | `${VCF_DIR}/` | Varies | Merged/concatenated VCF |
| `combined_for_pca.vcf.gz.csi` | `${VCF_DIR}/` | ~1-10 MB | CSI index (large genomes) |
| `combined_for_pca.stats.txt` | `${VCF_DIR}/` | ~10-100 KB | bcftools stats summary |

**stats.txt includes:**
- Total records, SNPs, indels
- Sample count
- Singleton counts
- Allele frequency distribution
- Ti/Tv ratios

---

## 4. Assessment: Room for Improvement?

### ✅ Strong Points

1. **Well-organized directory structure**
   - Clear separation by function (metrics, plots, PCA)
   - Configurable via environment variables
   - Dataset-specific SLURM log isolation

2. **Comprehensive SLURM logging**
   - Separate .output and .error files
   - Unique job identifiers in filenames
   - Fallback to /tmp if log directory unavailable
   - Dataset-specific subdirectories prevent collision

3. **Rich output documentation**
   - Multiple output formats (TSV, PNG, PDF)
   - Raw data + visualizations
   - Intermediate files preserved for debugging

4. **Flexible configuration**
   - All paths customizable
   - Sensible defaults
   - Environment variable overrides

### ⚠️ Areas for Improvement

#### 1. **SLURM Log Retention Strategy** (Minor)

**Issue:** No automatic cleanup or archival strategy for SLURM logs

**Current State:**
```bash
${LOG_BASE_PATH}/${dataset_name}/slurm/
  - 1D_12345678_0_job_12345678.output  # Never deleted
  - 1D_12345678_0_job_12345678.error
  - 1D_12345679_0_job_12345679.output  # Accumulates over time
  - 1D_12345679_0_job_12345679.error
  # ... potentially hundreds of log files ...
```

**Impact:** Low (logs are small, but can accumulate over months/years)

**Suggestions:**
- Add retention policy (e.g., "keep last 30 days")
- Archive old logs to compressed format
- Add cleanup script: `bash lib/cleanup_old_logs.sh --days 30`
- Document recommended cleanup schedule in README

#### 2. ~~**Granular Plot Regeneration**~~ ✅ FIXED (Feb 3, 2026)

~~**Issue:** All-or-nothing plot regeneration - missing 1 plot → regenerate all 17~~

**Status:** ✅ **IMPLEMENTED** - See `GRANULAR_PLOT_FIX.md`

Now regenerates only missing plots:
- ✅ Saves up to 8.5 minutes per re-run
- ✅ Reduced I/O and memory usage
- ✅ Improved logging shows which plots are regenerated
- ✅ Fully backward compatible

**Example:**
```bash
rm missingness_plots/Chr05*.png
bash step1d_interactive.sh --qc
# Before: Regenerate all 17 plots (~8.5 min)
# After:  Regenerate only Chr05 (~30 sec) ⚡
```

### 3. **Output Summary File** (Minor Enhancement)

**Issue:** No single "manifest" file listing all generated outputs

**Current State:** Users must check multiple directories to understand what was produced

**Suggestion:** Generate a summary file at workflow completion:

```bash
${WORK_DIR}/step1d_output_manifest.txt
---
Step1D Workflow Summary
Generated: 2026-02-03 14:23:45
Mode: PCA
Dataset: Apple_2024

Outputs:
✓ combined_for_pca.vcf.gz (1.2 GB, 1,234,567 variants, 150 samples)
✓ pca_analysis/pca.eigenvec (150 samples, 10 PCs)
✓ pca_analysis/pca_PC1_PC2.png
✓ pca_analysis/pca_scree.png
✓ king_duplicate_pairs.tsv (5 duplicate pairs detected)

SLURM Logs:
  /scratch/user/uqpha1/logs/Apple_2024/slurm/1D_12345678_*.{output,error}

Runtime: 2h 34m
```

**Benefits:**
- Quick verification of workflow completion
- Easy troubleshooting (missing files immediately visible)
- Useful for HPC job monitoring dashboards

#### 4. **SLURM Log Directory Naming** (Very Minor)

**Issue:** The subdirectory is always named `slurm/`, which can be generic if logs are stored in a shared location

**Current:**
```
/scratch/user/uqpha1/logs/Apple_2024/slurm/
/scratch/user/uqpha1/logs/Other_Dataset/slurm/
```

**Potential Improvement:**
```
/scratch/user/uqpha1/logs/Apple_2024/step1d/
/scratch/user/uqpha1/logs/Other_Dataset/step1d/
```

**Benefit:** More intuitive for multi-step pipelines (step1a, step1b, step1d, step2, etc.)

**Implementation:** Change `resolve_log_root()` call in `lib/slurm.sh` from `"slurm"` to `"step1d"`

#### 5. **Output Validation** (Enhancement)

**Issue:** No post-run validation that all expected outputs were created

**Current State:** Workflow completes successfully even if some plots fail to generate (non-fatal errors)

**Suggestion:** Add optional validation flag:
```bash
bash wrappers/sbatch/step1d_submit.sh Apple_2024 /data/vcfs --PCA --validate
```

**Validation checks:**
- All expected files exist
- Files are non-empty
- PNG/PDF files are valid images
- TSV files have expected column counts
- Exit with error if validation fails

#### 6. **Intermediate File Management** (Minor)

**Issue:** Some intermediate files are kept (e.g., PLINK2 .bed/.bim/.fam in PCA mode) without clear documentation

**Current State:**
- Useful for debugging
- Take up disk space
- Users may not know if they're safe to delete

**Suggestions:**
1. Add `--cleanup` flag to remove intermediate files after successful completion
2. Document which files are "final outputs" vs "intermediate" in README
3. Create separate `${STEP1D_PCA_DIR}/intermediate/` subdirectory

#### 7. **SLURM Log Filename Length** (Very Minor)

**Current Format:**
```
1D_12345678_0_Apple_GATK_1D_Apple_2024_12345678.output
```

**Issue:** Redundant job ID (appears twice: `%A` and `%j` are often the same)

**Suggestion:** Simplify to:
```
1D_Apple_2024_12345678.output
# Format: ${step_name}_${dataset_name}_%j.{output,error}
```

**Benefit:** Shorter filenames, easier to read at a glance

---

## Summary Assessment

### Overall: **Excellent** 🌟🌟

The current output management is **well-designed and production-ready**. After implementing granular plot regeneration, the system is even more efficient. Remaining suggestions are minor enhancements.

### Priority Recommendations

| Priority | Issue | Effort | Impact | Status |
|----------|-------|--------|--------|--------|
| ~~**High**~~ | ~~Granular plot regeneration~~ | ~~Low~~ | ~~High~~ | ✅ **FIXED** |
| **High** | Add output summary manifest | Low | High (improves UX) | Open |
| **Medium** | Document cleanup strategy | Low | Medium (prevents clutter) | Open |
| **Low** | Simplify SLURM log filenames | Low | Low (cosmetic) | Open |
| **Low** | Add `--validate` flag | Medium | Medium (catches errors) | Open |
| **Very Low** | Rename `slurm/` to `step1d/` | Low | Low (cosmetic) | Open |
| **Very Low** | Separate intermediate files | Low | Low (organization) | Open |

### Is It Perfect?

**Almost!** The current implementation is:
- ✅ Robust
- ✅ Well-documented
- ✅ Configurable
- ✅ Production-ready

The suggested improvements are **nice-to-haves** that would enhance usability and maintainability, but **none are critical**.

### Implementation Strategy

**Option 1: Leave as-is** (Recommended if pipeline is stable and working)
- Current design is solid
- Changes introduce testing overhead
- Focus on new features instead

**Option 2: Implement Remaining High-Priority Items** (Recommended if actively developing)
1. ~~Granular plot regeneration~~ ✅ **DONE** (Feb 3, 2026)
2. Add output manifest (30 minutes)
3. Document cleanup strategy in README (15 minutes)
4. Total remaining effort: ~45 minutes

**Option 3: Full Enhancement** (If time permits)
- Implement all 6 suggestions
- Total effort: ~4-6 hours
- Benefit: Slightly improved UX and maintainability

---

## Related Files

- **lib/slurm.sh**: SLURM log configuration (lines 82-104, 186-192)
- **lib/logging.sh**: Log directory resolution (lines 43-72)
- **config/environment.template.sh**: Default paths (lines 31, 28)
- **modules/step1d/README_FIRST.md**: User-facing documentation (lines 85-90)
- **modules/step1d/templates/master_vcf_analysis.sh**: Output generation logic
- **modules/step1d/templates/plink2_PCA.sh**: PCA-specific outputs
- **modules/step1d/GRANULAR_PLOT_FIX.md**: Details on the implemented plot regeneration fix

---

**Updates:**
- **Feb 3, 2026**: Implemented granular plot regeneration (Fix #2) - saves up to 8.5 minutes per re-run

---

**Conclusion:** Your Step1D module has excellent output management. After implementing granular plot regeneration, the system is even more efficient. Remaining suggestions are minor enhancements that would make a great system slightly better, but are not necessary for production use.
