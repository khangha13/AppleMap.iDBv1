# Step1D Function & Flag Map
**Report Date:** 2026-02-03  
**Purpose:** Comprehensive mapping of all entry points, flags, and overlaps

---

## Executive Summary

Step1D has **5 entry points** with **overlapping flags** and **inconsistent naming**. The system works but is confusing for users.

**Key Issues:**
- 🔴 Same flag names with different behaviors
- 🟡 Inconsistent capitalization (--PCA vs --pca vs --pca-only)
- 🟡 Hidden dependencies between scripts
- 🟢 Multiple ways to do the same thing (good flexibility, bad discoverability)

---

## 1. Entry Points Overview

| Entry Point | Type | Lines | Primary Use Case | User Facing |
|-------------|------|-------|------------------|-------------|
| **run_step1d.sh** | SLURM wrapper | 196 | Submit job to cluster | ✅ Yes |
| **step1d_interactive.sh** | Interactive wrapper | 514 | Local/interactive analysis | ✅ Yes |
| **master_vcf_analysis.sh** | Main analysis | 1,511 | Actual QC/PCA work | ⚠️ Advanced |
| **plink2_PCA.sh** | PCA helper | 507 | PCA computation | ❌ Internal |
| **prepare_combined_for_pca.sh** | VCF prep utility | 496 | Prepare VCF for PCA | ✅ Yes (new) |

---

## 2. Complete Flag Matrix

### 2.1 Comprehensive Flag Comparison

| Flag | run_step1d.sh | step1d_interactive.sh | master_vcf_analysis.sh | prepare_combined_for_pca.sh | plink2_PCA.sh |
|------|--------------|---------------------|----------------------|---------------------------|--------------|
| `--help`, `-h` | ✅ | ✅ | ✅ | ✅ | ❌ |
| `--dry-run`, `-n` | ✅ | ✅ | ✅ | ✅ | ❌ |
| `--force` | ❌ | ❌ | ❌ | ✅ | ❌ |
| `--beagle` | ✅ | ✅ | ✅ | ❌ | ❌ |
| `--qc` | ✅ | ✅ | ✅ | ❌ | ❌ |
| `--PCA` | ✅ | ⚠️ `--PCA` | ✅ `--PCA`/`--pca` | ❌ | ❌ |
| `--pca-only` | ❌ | ⚠️ Alias | ❌ | ❌ | ❌ |
| `--duplicate-check` | ✅ | ✅ | ✅ | ❌ | ❌ |
| `--remove-relatives` | ✅ | ✅ | ✅ | ❌ | ❌ |
| `--dir=PATH` | ❌ | ✅ | ❌ | ❌ (positional) | ❌ |
| `--vcf=LIST` | ❌ | ✅ | ❌ | ❌ | ❌ |

### 2.2 Flag Inconsistencies

#### Issue 1: Capitalization Chaos
```bash
# Different forms of the same flag:
run_step1d.sh --PCA              # Uppercase
run_step1d.sh --pca              # Lowercase (also accepted)
step1d_interactive.sh --PCA      # Uppercase
step1d_interactive.sh --pca      # Lowercase  
step1d_interactive.sh --pca-only # With dash (alias)
master_vcf_analysis.sh --PCA     # Uppercase
master_vcf_analysis.sh --pca     # Lowercase (also accepted)
```

**Status:** ⚠️ Confusing but functional (case-insensitive matching)

#### Issue 2: Mode Flags Only Work at Top Level
```bash
# These work:
run_step1d.sh /vcfs --PCA
step1d_interactive.sh --PCA

# These DON'T work (master_vcf_analysis.sh doesn't accept directory as arg):
master_vcf_analysis.sh /vcfs --PCA  # ❌ Wrong!

# Must use environment vars:
export VCF_DIR=/vcfs
bash master_vcf_analysis.sh --PCA   # ✅ Correct
```

**Status:** 🔴 Confusing - not documented clearly

---

## 3. Script-by-Script Deep Dive

### 3.1 run_step1d.sh (SLURM Wrapper)

**Purpose:** Submit Step1D analysis to SLURM cluster

**Signature:**
```bash
run_step1d.sh <dataset_name> <vcf_directory> [OPTIONS]
# OR
run_step1d.sh <vcf_directory> [OPTIONS]  # Dataset name = directory basename
```

**Flags:**
| Flag | Type | Default | Description | Conflicts |
|------|------|---------|-------------|-----------|
| `--beagle` | Boolean | false | Treat as Beagle-imputed VCFs | None |
| `--dry-run`, `-n` | Boolean | false | Preview without execution | None |
| `--qc` | Mode | ✅ (default) | Run QC metrics + plots | --PCA, --duplicate-check |
| `--PCA`, `--pca` | Mode | - | Run PCA analysis only | --qc, --duplicate-check |
| `--duplicate-check` | Mode | - | KING duplicate detection only | --qc, --PCA |
| `--remove-relatives` | Boolean | false | Remove relatives before PCA | Requires --PCA |

**Modes:** Mutually exclusive - only ONE of `--qc`, `--PCA`, or `--duplicate-check`

**What it does:**
1. Validates dataset/VCF directory
2. Reads Step1D config from pipeline_config.sh
3. Generates SLURM script from template
4. Submits job with `sbatch`
5. Returns job ID

**Example:**
```bash
run_step1d.sh NCBI /data/vcfs --PCA --remove-relatives
```

---

### 3.2 step1d_interactive.sh (Interactive Wrapper)

**Purpose:** Interactive/local analysis with user prompts

**Signature:**
```bash
step1d_interactive.sh [OPTIONS]
```

**Flags:**
| Flag | Type | Default | Description | Conflicts |
|------|------|---------|-------------|-----------|
| `--dir=PATH` | String | (prompts) | VCF directory path | None |
| `--vcf=LIST` | String | (auto-detect) | Comma-separated VCF prefixes | None |
| `--beagle` | Boolean | false | Beagle-imputed mode | None |
| `--qc` | Mode | ✅ (default) | QC-only mode | --PCA, --duplicate-check |
| `--PCA`, `--pca`, `--pca-only` | Mode | - | PCA-only mode | --qc, --duplicate-check |
| `--duplicate-check` | Mode | - | Duplicate detection only | --qc, --PCA |
| `--remove-relatives` | Boolean | false | Remove relatives before PCA | Requires --PCA |
| `--dry-run`, `-n` | Boolean | false | Preview only | None |
| `--help`, `-h` | Boolean | false | Show help | None |

**Special Features:**
- **VCF auto-detection:** Finds Chr*.vcf.gz files automatically
- **VCF filtering:** `--vcf=Chr01,Chr02` to analyze subset
- **Interactive prompts:** Guides user through selections
- **Fallback logic:** Tries *_snps.vcf.gz if Chr*.vcf.gz missing

**What it does:**
1. Prompts for VCF directory (unless `--dir=` provided)
2. Scans directory for VCF files
3. Validates chromosome count (expects 18 by default, 17 for Beagle)
4. Sets environment variables
5. Executes `master_vcf_analysis.sh` directly

**Example:**
```bash
# Interactive (prompts for directory):
step1d_interactive.sh --PCA

# Non-interactive:
step1d_interactive.sh --dir=/data/vcfs --vcf=Chr01,Chr02,Chr03 --qc
```

---

### 3.3 master_vcf_analysis.sh (Main Analysis Engine)

**Purpose:** Core analysis logic (QC metrics extraction + plotting + PCA)

**Signature:**
```bash
master_vcf_analysis.sh [OPTIONS]
```

**❗ Important:** Does NOT accept directory as argument! Uses environment variables.

**Flags:**
| Flag | Type | Default | Description | Conflicts |
|------|------|---------|-------------|-----------|
| `--dry-run`, `-n` | Boolean | false | Preview actions | None |
| `--beagle` | Boolean | false | Beagle-imputed metrics | None |
| `--qc` | Mode | ✅ (default) | QC metrics + plots | --PCA, --duplicate-check |
| `--PCA`, `--pca` | Mode | - | PCA analysis | --qc, --duplicate-check |
| `--duplicate-check` | Mode | - | Duplicate detection | --qc, --PCA |
| `--remove-relatives` | Boolean | false | Remove relatives (PCA) | Requires --PCA |
| `--help`, `-h` | Boolean | false | Show usage | None |

**Required Environment Variables:**
```bash
VCF_DIR=/path/to/vcfs              # REQUIRED
VCF_PATTERN="Chr%02d.vcf.gz"       # Optional, default shown
WORK_DIR=/path/to/output           # Optional, defaults to VCF_DIR
R_SCRIPTS_DIR=/path/to/Rscripts    # Optional, auto-detected
# ... many more (see config section in script)
```

**QC Mode Workflow:**
1. **Metrics Extraction** (per-chromosome)
   - Filter to biallelic SNPs: `bcftools view -m2 -M2 -v snps`
   - Extract site metrics (QUAL, QD, AC, AF, depth, etc.)
   - Generate per-chromosome TSV files
   - Combine into `variant_site_metrics.tsv`

2. **Plot Generation** (9 plot types):
   - Allele frequency distributions
   - Mean depth vs position
   - Missingness vs position
   - Depth vs missingness scatter
   - Site quality histograms
   - Heterozygosity rates
   - Quality-by-depth (QD)
   - Call rate heatmaps
   - Combined missingness grid

**PCA Mode Workflow:**
1. **VCF Preparation** (NEW - auto-triggers prepare_combined_for_pca.sh)
   - Detect merged VCF or concat chromosomes
   - Remove Chr00
   - Standardize to Chr01-Chr17
   - Generate bcftools stats

2. **PLINK2 Conversion**
   - Convert to `.pgen/.pvar/.psam` format
   - Filter: `--snps-only just-acgt --max-alleles 2`

3. **QC Filtering**
   - `--geno 0.05` (max 5% missing per variant)
   - `--mind 0.10` (max 10% missing per sample)
   - `--maf 0.01` (MAF ≥1%)

4. **Duplicate Detection** (if enabled)
   - KING kinship calculation
   - Flag pairs with kinship ≥0.45 (configurable)
   - Output: `king_duplicate_pairs.tsv`, `king_duplicate_samples.tsv`

5. **LD Pruning**
   - `--indep-pairwise 200 50 0.2` (200kb window, r²<0.2)

6. **PCA Computation**
   - 10 principal components
   - Variance explained per PC
   - Optional relative removal (KING cutoff 0.125)

7. **Plotting**
   - PC1 vs PC2 scatter
   - Duplicate samples flagged in red
   - Variance explained in axis labels

**Duplicate-Check Mode:**
- Same as PCA steps 1-4, then stops
- Outputs only the duplicate tables

**Beagle Mode Differences:**
- Different metrics: AF, DR2 (dosage R²), IMP (imputation flag)
- Skips depth-based metrics (MEAN_DEPTH, QD)
- Skips depth-based plots

---

### 3.4 plink2_PCA.sh (PCA Helper - Internal)

**Purpose:** Execute PLINK2 PCA workflow (called by master_vcf_analysis.sh)

**Signature:**
```bash
plink2_PCA.sh <vcf_dir> <rscripts_dir> [plink2_bin] [bcftools_bin] \
              [remove_relatives] [show_labels] [label_size] [use_ggrepel] \
              [merged_vcf_pattern] [duplicate_mode] [duplicate_king_threshold] \
              [run_mode]
```

**12 positional parameters!** 😱

**Parameters:**
| Position | Name | Type | Default | Description |
|----------|------|------|---------|-------------|
| 1 | vcf_dir | String | REQUIRED | Source VCF directory |
| 2 | rscripts_dir | String | REQUIRED | R scripts location |
| 3 | plink2_bin | String | `plink2` | PLINK2 binary path |
| 4 | bcftools_bin | String | `bcftools` | bcftools binary path |
| 5 | remove_relatives | Boolean | `false` | Remove relatives (KING) |
| 6 | show_labels | Boolean | `true` | Show sample IDs on plot |
| 7 | label_size | Number | `3` | Plot label size |
| 8 | use_ggrepel | Boolean | `true` | Use ggrepel for labels |
| 9 | merged_vcf_pattern | String | `*merged*.vcf.gz` | Merged VCF glob |
| 10 | duplicate_mode | String | `flag` | `off`/`flag`/`remove` |
| 11 | duplicate_king_threshold | Number | `0.45` | Kinship threshold |
| 12 | run_mode | String | `pca` | `pca` or `duplicate` |

**No flag interface - pure positional!**

**What it does:**
1. Resolve VCF input (merged vs per-chromosome)
2. Convert to PLINK2 format
3. Apply QC filters
4. Run duplicate detection (if enabled)
5. LD pruning
6. PCA computation
7. Generate plots
8. Cleanup intermediate files

**Status:** 🟡 Internal script, but complex interface

---

### 3.5 prepare_combined_for_pca.sh (VCF Prep Utility - NEW)

**Purpose:** Prepare standardized VCF for PCA input

**Signature:**
```bash
prepare_combined_for_pca.sh <vcf_directory> [OPTIONS]
```

**Flags:**
| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--force` | Boolean | false | Overwrite existing output |
| `--dry-run` | Boolean | false | Preview without creating files |
| `--help`, `-h` | Boolean | false | Show help message |

**Output:**
- `combined_for_pca.vcf.gz` - Combined, cleaned VCF
- `combined_for_pca.vcf.gz.csi` - Index
- `combined_for_pca.stats.txt` - bcftools stats (SNP counts!)

**What it does:**
1. **VCF Selection Logic:**
   - Prefer merged VCF: `*merged*.vcf.gz` or `*merge*.vcf.gz` (exclude names with 'Chr')
   - Fallback: Concatenate all `*.vcf.gz` (exclude `Chr00*.vcf.gz`)

2. **Chr00 Removal:**
   - If merged VCF: `bcftools view -t ^Chr00`
   - If concatenating: Skip `Chr00*.vcf.gz` files

3. **Chromosome Standardization:**
   - Map: `1→Chr01`, `chr1→Chr01`, `Chr1→Chr01`, etc.
   - Range: Chr01-Chr17
   - Use: `bcftools annotate --rename-chrs`

4. **Tabix Workaround:**
   - Run `tabix -f -p vcf` before renaming (per past problems §contig undefined)

5. **Sort & Index:**
   - `bcftools sort -Oz`
   - `bcftools index -f -c` (CSI format)

6. **Generate Stats:**
   - `bcftools stats > combined_for_pca.stats.txt`
   - Extract SNP count in summary

**Example:**
```bash
prepare_combined_for_pca.sh /data/vcfs
# Creates:
#   /data/vcfs/combined_for_pca.vcf.gz
#   /data/vcfs/combined_for_pca.vcf.gz.csi
#   /data/vcfs/combined_for_pca.stats.txt
```

**Integration:**
- NEW: `master_vcf_analysis.sh --PCA` now auto-calls this if combined VCF missing/outdated
- Standalone: Can be run independently

---

## 4. Overlap Analysis

### 4.1 Duplicate Functionality

#### Overlap 1: VCF Combination Logic

**Where duplicated:**
- `plink2_PCA.sh` (lines 167-267)
- `prepare_combined_for_pca.sh` (entire script)

**What overlaps:**
- Merged VCF detection pattern matching
- Fallback to per-chromosome concatenation
- VCF file discovery logic

**Why it exists:**
- `plink2_PCA.sh` originally did this inline
- `prepare_combined_for_pca.sh` extracted it to standalone utility
- `plink2_PCA.sh` still has old logic (not yet removed)

**Impact:** 🟡 Moderate - could diverge if one is updated

**Resolution:** Make `plink2_PCA.sh` always expect pre-combined VCF

---

#### Overlap 2: Dry-Run Mode

**Where duplicated:**
- `run_step1d.sh` (passes to master)
- `step1d_interactive.sh` (passes to master)
- `master_vcf_analysis.sh` (implements it)
- `prepare_combined_for_pca.sh` (implements it)

**What overlaps:**
- Flag parsing: `--dry-run` or `-n`
- Preview logic without file creation

**Why it exists:**
- Each entry point needs to pass it down
- Each worker script needs to implement it

**Impact:** ✅ Low - consistent behavior

---

#### Overlap 3: Mode Selection (QC/PCA/Duplicate)

**Where duplicated:**
- `run_step1d.sh` (parses, passes to master)
- `step1d_interactive.sh` (parses, passes to master)
- `master_vcf_analysis.sh` (implements modes)
- `plink2_PCA.sh` (has `run_mode` parameter)

**What overlaps:**
- Mutual exclusivity checking
- Mode validation
- Error messages

**Impact:** 🟡 Moderate - error messages differ slightly

---

### 4.2 Missing Functionality Overlap

**Gap 1: No unified config reader**
- Each script parses environment variables independently
- No shared validation of VCF_DIR, VCF_PATTERN, etc.

**Gap 2: No shared VCF discovery**
- `step1d_interactive.sh` has complex logic to find VCFs
- `master_vcf_analysis.sh` has different logic
- `prepare_combined_for_pca.sh` has third implementation

**Gap 3: No shared dry-run implementation**
- Each script reimplements "preview without execution"
- Inconsistent verbosity levels

---

## 5. Flag Consistency Issues

### 5.1 Critical Issues

#### Issue 1: Capitalization Inconsistency

| Flag | Variations Found | Standard? |
|------|------------------|-----------|
| PCA | `--PCA`, `--pca`, `--pca-only` | ❌ No |
| Dry-run | `--dry-run`, `-n` | ✅ Yes |
| Help | `--help`, `-h` | ✅ Yes |
| Beagle | `--beagle` | ✅ Yes |

**Recommendation:** Standardize on lowercase with case-insensitive matching:
- `--pca` (canonical)
- `--PCA` (accepted as alias)

---

#### Issue 2: Positional vs. Flag Arguments

| Script | Positional Args | Flag Args | Style |
|--------|----------------|-----------|-------|
| run_step1d.sh | 2 (dataset, vcf_dir) | 6 flags | ✅ Mixed |
| step1d_interactive.sh | 0 | 9 flags | ✅ Flags only |
| master_vcf_analysis.sh | 0 | 7 flags | ✅ Flags only |
| prepare_combined_for_pca.sh | 1 (vcf_dir) | 3 flags | ✅ Mixed |
| plink2_PCA.sh | **12** | 0 | ❌ Positional only! |

**Recommendation:** Convert `plink2_PCA.sh` to flag-based interface

---

#### Issue 3: Mode Flags Not Orthogonal

**Current behavior:**
```bash
# These are mutually exclusive:
--qc
--PCA  
--duplicate-check

# This is a modifier (only works with --PCA):
--remove-relatives
```

**Confusing aspect:** `--remove-relatives` is a flag, but depends on mode

**Recommendation:** Make dependency explicit in help text (already done) ✅

---

### 5.2 Minor Inconsistencies

#### Inconsistency 1: Help Flag Naming

```bash
run_step1d.sh        # No --help (only usage message on error)
step1d_interactive.sh --help, -h  ✅
master_vcf_analysis.sh --help, -h  ✅
prepare_combined_for_pca.sh --help, -h  ✅
plink2_PCA.sh        # No --help ❌
```

#### Inconsistency 2: Boolean Normalization

**Different patterns used:**
- `normalize_bool()` in master_vcf_analysis.sh
- `normalize_bool()` in plink2_PCA.sh (duplicated code!)
- Direct boolean checks in prepare_combined_for_pca.sh

---

## 6. Configuration Hell

### 6.1 Environment Variables (master_vcf_analysis.sh)

**28+ environment variables!**

**Core inputs:**
```bash
VCF_DIR          # Input directory (REQUIRED)
VCF_PATTERN      # File naming pattern
WORK_DIR         # Output directory
R_SCRIPTS_DIR    # R script location
```

**Output directories (all optional):**
```bash
MEAN_DP_TSV               # Default: SNP_site_meanDP.tsv
SITE_METRICS_TSV          # Default: variant_site_metrics.tsv
METRICS_BY_CHROM_DIR      # Default: site_metrics_per_chromosome
DEPTH_PLOTS_DIR           # Default: depth_plots
MISSINGNESS_PLOTS_DIR     # Default: missingness_plots
DEPTH_MISSINGNESS_DIR     # Default: depth_vs_missingness
SITE_QUALITY_DIR          # Default: site_quality_plots
HETEROZYGOSITY_DIR        # Default: heterozygosity_plots
QD_PLOTS_DIR              # Default: quality_by_depth_plots
CALL_RATE_HEATMAP_DIR     # Default: call_rate_heatmaps
AF_PLOTS_DIR              # Default: af_distribution_plots
```

**Plot settings:**
```bash
PLOT_IMAGE_FORMAT         # Default: png
CALL_RATE_HEATMAP_BINS    # Default: 100
AF_HIST_BINS              # Default: 50
```

**PCA settings (from pipeline_config.sh):**
```bash
STEP1D_PCA_DIR                  # Default: pca_analysis
STEP1D_REMOVE_RELATIVES         # Default: false
STEP1D_PCA_SHOW_LABELS          # Default: true
STEP1D_PCA_LABEL_SIZE           # Default: 1.5
STEP1D_PCA_USE_GGREPEL          # Default: true
STEP1D_PCA_MERGED_PATTERN       # Default: *merged*.vcf.gz,*merge*.vcf.gz
STEP1D_PCA_FORCE_CONCAT         # Default: false
STEP1D_PCA_MERGED_EXCLUDE_CHR   # Default: true
STEP1D_DUPLICATE_MODE           # Default: flag
STEP1D_DUPLICATE_KING_THRESHOLD # Default: 0.45
```

**Status:** 🔴 Overwhelming - needs documentation

---

## 7. Recommendations

### 7.1 High Priority (Usability)

#### 1. Unify Flag Naming
**Problem:** `--PCA`, `--pca`, `--pca-only` all mean the same thing

**Solution:**
```bash
# Standardize on lowercase, accept uppercase as alias
--pca          # Canonical
--PCA          # Accepted (case-insensitive)
--pca-only     # Remove (confusing)
```

**Effort:** 15 minutes  
**Impact:** High (reduces confusion)

---

#### 2. Create Quick Reference Guide
**Problem:** Users don't know which entry point to use

**Solution:** Create `QUICK_START.md`:
```markdown
# Step1D Quick Start

## I want to...

### Run QC on chromosome VCFs (SLURM):
bash modules/step1d/bin/run_step1d.sh /path/to/vcfs

### Run PCA on chromosome VCFs (SLURM):
bash modules/step1d/bin/run_step1d.sh /path/to/vcfs --PCA

### Run QC interactively:
bash wrappers/interactive/step1d_interactive.sh

### Just prepare a combined VCF:
bash modules/step1d/bin/prepare_combined_for_pca.sh /path/to/vcfs
```

**Effort:** 30 minutes  
**Impact:** High (first-time user experience)

---

#### 3. Add Helpful Error Messages
**Problem:** Wrong usage produces cryptic errors

**Solution:** Detect common mistakes:
```bash
# If user tries:
bash master_vcf_analysis.sh /path/to/vcfs --PCA

# Show:
"❌ master_vcf_analysis.sh doesn't accept directory arguments.
Use: export VCF_DIR=/path/to/vcfs
     bash master_vcf_analysis.sh --PCA

Or use a wrapper:
     bash bin/run_step1d.sh /path/to/vcfs --PCA"
```

**Effort:** 1 hour  
**Impact:** Medium (prevents frustration)

---

### 7.2 Medium Priority (Maintainability)

#### 4. Refactor plink2_PCA.sh Interface
**Problem:** 12 positional parameters is unmaintainable

**Solution:** Convert to flags:
```bash
# Before:
plink2_PCA.sh $vcf_dir $rscripts_dir $plink2 $bcftools $remove_rel...

# After:
plink2_PCA.sh --vcf-dir=/path --rscripts-dir=/path \
              --remove-relatives --duplicate-mode=flag
```

**Effort:** 2-3 hours  
**Impact:** Medium (internal script)

---

#### 5. Consolidate VCF Discovery Logic
**Problem:** 3 different implementations of "find VCF files"

**Solution:** Create `lib/vcf_discovery.sh`:
```bash
find_vcfs() {
    local dir="$1"
    local pattern="${2:-*.vcf.gz}"
    # Single implementation used everywhere
}
```

**Effort:** 2 hours  
**Impact:** Medium (reduces duplication)

---

### 7.3 Low Priority (Nice to Have)

#### 6. Create Unified Config Validator
**Problem:** Each script validates config independently

**Solution:** `lib/step1d_config.sh`:
```bash
validate_step1d_config() {
    # Check VCF_DIR exists
    # Check R_SCRIPTS_DIR exists
    # Validate all env vars
    # Return meaningful errors
}
```

**Effort:** 2 hours  
**Impact:** Low (preventive)

---

#### 7. Environment Variable Documentation
**Problem:** 28+ env vars with no central docs

**Solution:** Generate `CONFIG_REFERENCE.md` from script:
```bash
# Extract all ${VAR:-default} patterns
# Auto-generate docs table
```

**Effort:** 1 hour  
**Impact:** Low (reference material)

---

## 8. Summary Table: What to Use When

| Scenario | Use This | Command |
|----------|----------|---------|
| **QC on cluster** | run_step1d.sh | `run_step1d.sh /vcfs` |
| **PCA on cluster** | run_step1d.sh | `run_step1d.sh /vcfs --PCA` |
| **Interactive QC** | step1d_interactive.sh | `step1d_interactive.sh` |
| **Prepare VCF only** | prepare_combined_for_pca.sh | `prepare_combined_for_pca.sh /vcfs` |
| **Advanced/custom** | master_vcf_analysis.sh | `export VCF_DIR=/vcfs; master_vcf_analysis.sh --PCA` |
| **NEVER use directly** | plink2_PCA.sh | ❌ Internal only |

---

## 9. Overlap Summary

### Intentional Duplication (OK)
- ✅ Dry-run flag parsing (each entry point needs it)
- ✅ Mode validation (wrappers validate, worker implements)
- ✅ Help text (each script has its own)

### Unintentional Duplication (Fix It)
- 🔴 VCF combination logic (plink2_PCA.sh vs prepare_combined_for_pca.sh)
- 🔴 Boolean normalization function (copy-pasted)
- 🟡 VCF discovery patterns (3 implementations)

### Missing Shared Code (Add It)
- Config validation
- VCF discovery helper
- Common error messages

---

## 10. User Experience Rating

| Entry Point | Discoverability | Ease of Use | Documentation | Overall |
|-------------|----------------|-------------|---------------|---------|
| run_step1d.sh | ⭐⭐⭐ | ⭐⭐⭐⭐ | ⭐⭐⭐ | **⭐⭐⭐** |
| step1d_interactive.sh | ⭐⭐⭐⭐ | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐ | **⭐⭐⭐⭐** |
| master_vcf_analysis.sh | ⭐⭐ | ⭐⭐⭐ | ⭐⭐⭐ | **⭐⭐** |
| prepare_combined_for_pca.sh | ⭐⭐⭐ | ⭐⭐⭐⭐ | ⭐⭐⭐⭐ | **⭐⭐⭐** |
| plink2_PCA.sh | ⭐ | ⭐ | ⭐ | **⭐** |

**Best for users:** `step1d_interactive.sh` (5/5 ease of use)  
**Most powerful:** `master_vcf_analysis.sh` (but hardest to use)  
**Needs improvement:** `plink2_PCA.sh` (internal, but still...)

---

## Appendix A: Complete Flag Reference

See full table in §2.1

## Appendix B: Environment Variable Reference

(Would generate from script - too long for this report)

## Appendix C: Related Documentation

- `README_FIRST.md` - Module overview
- `HOW_TO_RUN_PCA.md` - PCA-specific guide
- `DIRECT_CALL_GUIDE.md` - Advanced usage
- `VCF_Analysis_Quick_Reference.md` - Command cheat sheet
- `FILTERING_OPTIONS_EXPLAINED.md` - Parameter details

---

**Report compiled by:** AI Assistant  
**Next steps:** Review recommendations §7 and prioritize fixes
