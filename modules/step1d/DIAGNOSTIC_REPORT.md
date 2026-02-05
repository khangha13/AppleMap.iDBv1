# Step1D Module - Comprehensive Diagnostic Report
**Generated:** 2026-02-03  
**Module Version:** 1.0  
**Total Lines of Code:** ~3,497 (bash + R)

---

## Executive Summary

Step1D is a **well-architected VCF analysis module** that combines QC metrics extraction, visualization, PCA analysis, and duplicate detection into a unified workflow. The module demonstrates **strong engineering practices** with robust error handling, extensive configuration options, and comprehensive documentation. However, there are opportunities for improvement in consistency, code duplication, and testing infrastructure.

**Overall Rating:** ⭐⭐⭐⭐ (4/5)

---

## 1. Architecture & Design

### ✅ Strong Points

#### 1.1 Modular Design
- **Clean separation of concerns:**
  - `bin/` - Entry points and utilities
  - `templates/` - Core workflow scripts
  - `Rscripts/` - Visualization layer
  - Documentation files at module root
- **Single Responsibility:** Each script has a well-defined purpose
- **Composability:** Scripts can be used independently or orchestrated together

#### 1.2 Configuration Management
- **Hierarchical configuration system:**
  ```
  Pipeline Config → Environment Variables → Script Defaults
  ```
- **20+ configurable parameters** for Step1D (all prefixed with `STEP1D_`)
- **Centralized defaults** in `config/pipeline_config.sh`
- **Override flexibility** at multiple levels (environment, CLI flags)

#### 1.3 Integration with Pipeline Infrastructure
- **Leverages shared libraries:**
  - `lib/logging.sh` - Standardized logging
  - `lib/slurm.sh` - SLURM job management
  - `config/pipeline_config.sh` - Configuration system
- **PIPELINE_ROOT resolution** handles both interactive and SLURM execution contexts
- **Graceful fallback** when module system unavailable

#### 1.4 SLURM Integration
- **Proper job submission** via `run_step1d.sh` wrapper
- **Resource configuration** specific to Step1D workloads (6 CPUs, 80GB RAM, 72h)
- **Dynamic script generation** from templates
- **Job metadata tracking** (dataset name, job ID)

### ⚠️ Weaknesses

#### 1.5 Code Duplication
- **PIPELINE_ROOT resolution logic** appears in 4+ files with slight variations:
  ```bash
  # Pattern repeated across:
  # - run_step1d.sh
  # - master_vcf_analysis.sh
  # - plink2_PCA.sh
  # - prepare_combined_for_pca.sh
  ```
  **Impact:** Maintenance burden, potential for divergence
  **Recommendation:** Extract to `lib/path_resolution.sh`

#### 1.6 Inconsistent Error Handling
- **master_vcf_analysis.sh:** Uses `trap` with logging + `set -Eeuo pipefail`
- **plink2_PCA.sh:** Uses `trap` with logging + `set -Eeuo pipefail`
- **prepare_combined_for_pca.sh:** Uses `set -euo pipefail` (no trap, no `-E`)
- **run_step1d.sh:** No explicit error handling flags

**Impact:** Inconsistent error propagation behavior  
**Recommendation:** Standardize on `set -Eeuo pipefail` + trap pattern everywhere

#### 1.7 Module Loading Strategy
- **Hard-coded module names** in `master_vcf_analysis.sh`:
  ```bash
  module load miniforge/25.3.0-3
  module load bcftools/1.18-gcc-12.3.0
  module load plink/2.00a3.6-gcc-11.3.0
  ```
- **Inconsistent approach:** `prepare_combined_for_pca.sh` uses generic `module load bcftools`
- **No configuration variable** for conda environment name in main script

**Impact:** Brittle across HPC environments, version lock-in  
**Recommendation:** Move module specifications to `pipeline_config.sh`

---

## 2. Code Quality & Standards

### ✅ Strong Points

#### 2.1 Bash Best Practices
- **Comprehensive safety flags:** Most scripts use `set -Eeuo pipefail`
- **Path handling:** Proper use of `"${VAR}"` quoting
- **Array usage:** Modern bash arrays (`declare -a`, `mapfile`)
- **Subshell isolation:** PCA helper runs in `(cd "${PCA_OUTPUT_DIR}" && ...)`
- **Temporary files:** Proper use of `mktemp` with TMPDIR

#### 2.2 R Code Quality
- **Consistent style:** All R scripts follow similar patterns
- **Package management:**
  - Graceful checks: `requireNamespace(..., quietly = TRUE)`
  - Clear error messages when packages missing
- **Common utilities:** Shared `common_plot_utils.R` reduces duplication
- **Data validation:** Input files and columns verified before processing
- **Modern libraries:** Uses `data.table` (fast), `ggplot2` (publication-quality), `ragg` (better rendering)

#### 2.3 Documentation
- **Comprehensive user guides:**
  - `README_FIRST.md` - Overview and quick start
  - `HOW_TO_RUN_PCA.md` - Detailed PCA workflow
  - `DIRECT_CALL_GUIDE.md` - Interactive usage
  - `FILTERING_OPTIONS_EXPLAINED.md` - Parameter reference
  - `VCF_Analysis_Quick_Reference.md` - Command cheat sheet
- **Inline documentation:** Header blocks in all scripts
- **Usage functions:** Help text available via `--help` flag

#### 2.4 Logging
- **Structured logging levels:** DEBUG, INFO, WARN, ERROR, FATAL
- **Colored output:** User-friendly terminal display
- **Dual logging:** Console + file (when in pipeline mode)
- **Context preservation:** Dataset name, module name included in logs
- **Stderr routing:** Prevents log pollution in command substitution

### ⚠️ Weaknesses

#### 2.5 Variable Naming Inconsistencies
```bash
# Inconsistent prefix patterns:
STEP1D_PCA_ONLY          # STEP1D_ prefix
PCA_GENO                 # No prefix
REMOVE_RELATIVES_CLI     # _CLI suffix
REMOVE_RELATIVES_RAW     # _RAW suffix
REMOVE_RELATIVES         # No suffix
```

**Impact:** Harder to trace variable origins  
**Recommendation:** Standardize on:
- `STEP1D_*` for configuration
- `*_CLI` for command-line arguments
- `*_RAW` for pre-normalized values
- No suffix for normalized/final values

#### 2.6 Magic Numbers & Strings
- **Hard-coded thresholds** scattered throughout:
  ```bash
  --geno 0.05    # Should be configurable
  --mind 0.10
  --maf 0.01
  --king-cutoff 0.125
  --indep-pairwise 200 50 0.2
  ```
- **Some configurable, others not:**
  - `STEP1D_DUPLICATE_KING_THRESHOLD` - ✅ Configurable
  - LD pruning parameters - ❌ Hard-coded

**Recommendation:** Make all QC thresholds configurable via environment variables

#### 2.7 Incomplete Input Validation
- **master_vcf_analysis.sh:** Validates VCF existence, checks commands
- **plink2_PCA.sh:** Validates directories, checks tools
- **prepare_combined_for_pca.sh:** Basic directory checks
- **Missing validations:**
  - Numeric parameter ranges (e.g., `PCA_GENO` should be 0-1)
  - File format verification (VCF headers)
  - Disk space checks before large operations

---

## 3. Feature Completeness

### ✅ Strong Points

#### 3.1 Comprehensive QC Metrics
- **Site-level metrics extraction:**
  - Quality scores (QUAL, QD)
  - Allele frequencies (AC, AF)
  - Depth statistics (MEAN_DEPTH, DP_NON_MISSING)
  - Missingness (CALL_RATE, MISSING_RATE)
  - Hardy-Weinberg (InbreedingCoeff, ExcessHet)
  - Heterozygosity rates
- **Beagle mode support:** Alternative metrics (AF, DR2, IMP) for imputed data

#### 3.2 Rich Visualization Suite
- **9 plot types:**
  1. Allele frequency distributions (per-chr + combined)
  2. Mean depth vs position
  3. Missingness vs position
  4. Depth vs missingness scatter
  5. Site quality histograms
  6. Heterozygosity rates
  7. Quality-by-depth (QD)
  8. Call rate heatmaps
  9. PCA scatter plots (with duplicate flagging)
- **Publication-ready:** High DPI (300), vector-friendly formats
- **Customizable:** Image format, bins, label sizes configurable

#### 3.3 PCA Workflow
- **Complete pipeline:**
  1. VCF combination (merged or concatenated)
  2. Conversion to PLINK2 format
  3. QC filtering (geno, mind, maf)
  4. LD pruning
  5. Optional duplicate removal
  6. Optional relative removal
  7. PCA computation (10 PCs + variant weights)
  8. Visualization with variance explained
- **Duplicate detection:** KING kinship-based (configurable: off/flag/remove)
- **Relative removal:** Configurable KING cutoff

#### 3.4 **NEW: VCF Preparation Utility** ✨
- **`prepare_combined_for_pca.sh`** - Recently integrated
  - Auto-detects merged VCF or concatenates per-chromosome files
  - Chromosome name standardization (→ `ChrNN`)
  - `Chr00` filtering
  - **bcftools stats generation** (includes SNP counts) 
  - Automatic integration with `--PCA` workflow
  - Timestamp-based regeneration logic

**Impact:** Major improvement! Reduces manual preprocessing and ensures consistent input format.

### ⚠️ Weaknesses

#### 3.5 Missing Features

##### Sample-Level QC
- **Current:** Site-level metrics only
- **Missing:** Per-sample statistics (missing rate, heterozygosity, relatedness)
- **Impact:** Cannot identify problematic samples before analysis
- **Recommendation:** Add `--sample-qc` mode

##### Variant Filtering Options
- **Current:** Hard-coded biallelic SNPs filter (`-m2 -M2 -v snps`)
- **Missing:** User-configurable filters (e.g., include indels, multiallelic sites)
- **Impact:** Limits flexibility for different analysis types
- **Recommendation:** Add `STEP1D_VARIANT_FILTER` configuration

##### Sex Chromosome Handling
- **Current:** Assumes autosomes (Chr01-Chr17)
- **Missing:** Optional sex chromosome inclusion
- **Impact:** Cannot analyze full genome
- **Note:** May be intentional for this species/project

##### Resume/Checkpoint System
- **Current:** Some reuse logic (e.g., existing `*_snps.vcf.gz`)
- **Missing:** Comprehensive checkpointing
- **Impact:** Restarting after failure repeats expensive steps
- **Recommendation:** Add `--resume` flag with state tracking

---

## 4. Consistency Analysis

### 🔍 Cross-File Consistency

#### 4.1 Naming Conventions

| Aspect | Status | Details |
|--------|--------|---------|
| **File naming** | ✅ Consistent | `snake_case.sh`, `PascalCase.R` |
| **Function naming** | ✅ Consistent | `snake_case()` in bash, `camelCase()` in R |
| **Variable naming** | ⚠️ Mixed | See 2.5 above |
| **Constants** | ✅ Consistent | `UPPER_SNAKE_CASE` |

#### 4.2 Configuration Patterns

**Consistent:**
- All Step1D configs use `STEP1D_` prefix
- Dual-source pattern: `${ENV_VAR:-default}` everywhere
- Boolean normalization: `normalize_bool()` function reused

**Inconsistent:**
- Some scripts use `_RAW` suffix, others don't
- Numeric normalization done inline vs. function (`normalize_num()` only in plink2_PCA)

#### 4.3 Error Handling Patterns

**Consistent:**
- All scripts use `log_error()` before `exit 1`
- File existence checks use `check_file()` or `[ ! -f ]`

**Inconsistent:**
- Trap usage (see 1.6)
- Error message format varies:
  ```bash
  log_error "Message"           # Most common
  echo "❌ Message" >&2; exit 1  # In arg parsing
  ```

#### 4.4 Dependency Checking

**Consistent:**
- R scripts all check packages the same way
- Bash scripts use `command -v` for tool checks

**Inconsistent:**
- Module loading done differently per script
- Some scripts attempt `module load`, others assume tools in PATH

---

## 5. Integration Quality

### ✅ Strong Points

#### 5.1 With Pipeline Infrastructure
- **Seamless config inheritance** from `pipeline_config.sh`
- **Logging integration** via `lib/logging.sh`
- **SLURM job management** via `lib/slurm.sh`
- **Path resolution** handles SLURM spool directory issues

#### 5.2 With External Tools
- **bcftools:**
  - Extensive use (view, query, concat, annotate, index, stats)
  - Graceful handling of missing formats (VCF vs. BCF)
  - Proper indexing (TBI for VCF, CSI for large chromosomes)
- **PLINK2:**
  - Modern `--pfile` format
  - Comprehensive QC options
  - KING kinship computation
- **R ecosystem:**
  - Modern data.table (fast I/O)
  - ggplot2 (publication-quality plots)
  - ragg (better font rendering)

#### 5.3 Data Flow
```
VCF Files → Site Metrics TSV → R Plots
          ↓
          → Combined VCF → PLINK2 → PCA → R Plots
```
**Clean pipeline:** Each stage produces reusable intermediate files

### ⚠️ Weaknesses

#### 5.4 Tight Coupling to HPC Environment
- **Hard-coded module names** (see 1.7)
- **Assumes SLURM** (though can run interactively)
- **TMPDIR dependency** (fails gracefully, but assumption present)

**Impact:** Portability issues on non-HPC systems  
**Recommendation:** Add Docker/Singularity container option

#### 5.5 No Version Checking
- **Tools:** No verification of bcftools/plink2/R versions
- **Modules:** Loads specific versions but doesn't verify compatibility
- **Packages:** R packages checked for presence, not version

**Impact:** Silent failures with incompatible versions  
**Recommendation:** Add version checks in initialization

---

## 6. Testing & Validation

### ⚠️ Major Gap

#### 6.1 No Automated Tests
- **No test suite** present
- **No fixtures** for validation
- **No CI/CD** integration

**Impact:** Regression risk, manual validation burden

#### 6.2 Testing Recommendations

**Unit Tests:**
```bash
# test/step1d/test_prepare_combined.sh
# - Test Chr00 removal
# - Test chromosome renaming
# - Test dry-run mode
# - Test with missing files
```

**Integration Tests:**
```bash
# test/step1d/test_full_workflow.sh
# - QC-only mode
# - PCA-only mode
# - Duplicate-check mode
# - Beagle mode
```

**Fixtures:**
```
test/fixtures/step1d/
  ├── minimal_vcf/          # 2 samples, 100 SNPs
  ├── beagle_vcf/           # Imputed data
  └── expected_outputs/     # Reference outputs
```

---

## 7. Performance Characteristics

### ✅ Strong Points

#### 7.1 Parallelization Opportunities
- **Per-chromosome operations** (metrics extraction) naturally parallel
- **Independent plot generation** can be parallelized
- **PLINK2** inherently multi-threaded (6 CPUs allocated)

#### 7.2 Memory Efficiency
- **Streaming operations:** bcftools processes VCFs without loading into memory
- **R data.table:** Fast, memory-efficient data loading
- **Temp file cleanup:** Intermediate files removed after use

#### 7.3 Incremental Processing
- **Reuses existing files:**
  - `*_snps.vcf.gz` (filtered SNPs)
  - Per-chromosome metrics TSV
  - Combined VCF for PCA
  - PLINK2 intermediate files (with param-based invalidation)
- **Smart regeneration:** Only recreates when:
  - Source files newer
  - Parameters changed
  - Output missing

### ⚠️ Weaknesses

#### 7.4 No Explicit Parallelization
- **Sequential processing:** Per-chromosome metrics extracted one at a time
- **Opportunity:** Could use `parallel` or `xargs -P` for per-chromosome work

#### 7.5 Large Memory Footprint
- **80GB RAM allocation** - generous but necessary for large VCFs
- **R plotting:** Loads entire metrics TSV into memory
- **Potential issue:** Very large datasets (>100M variants) may hit memory limits

---

## 8. Security & Robustness

### ✅ Strong Points

#### 8.1 Input Sanitization
- **Path validation:** Checks directory/file existence before use
- **Numeric validation:** `normalize_num()` ensures valid numbers
- **Boolean normalization:** Prevents shell injection via boolean flags

#### 8.2 Temp File Handling
- **Unique temp directories:** `mktemp -d` with random suffix
- **Cleanup on success:** Removes temp files explicitly
- **TMPDIR usage:** Respects HPC scratch space

### ⚠️ Weaknesses

#### 8.3 No Cleanup on Failure
- **Trap handling incomplete:** Some scripts don't clean up temp files on error
- **Partial outputs:** Failed runs may leave broken files

**Recommendation:**
```bash
cleanup() {
  rm -rf "${TEMP_DIR}" 2>/dev/null || true
}
trap cleanup EXIT ERR
```

#### 8.4 No Input Size Limits
- **Unbounded operations:** Will attempt to process arbitrarily large files
- **No disk space checks** before concatenation
- **Risk:** Out-of-space errors mid-operation

---

## 9. Maintainability

### ✅ Strong Points

#### 9.1 Self-Documenting Code
- **Descriptive variable names:** `COMBINED_VCF`, `PCA_OUTPUT_DIR`
- **Inline comments:** Explain "why" not just "what"
- **Section markers:** Clear visual separation in long scripts

#### 9.2 Version Control Friendly
- **No generated files** in repo
- **Template-based approach:** Core logic in version control
- **Configuration externalized:** Environment-specific values separate

#### 9.3 Extensibility
- **Modular R scripts:** Easy to add new plot types
- **Hook points:** Clear places to add preprocessing steps
- **Configuration-driven:** New features can be feature-flagged

### ⚠️ Weaknesses

#### 9.4 Long Functions
- **`master_vcf_analysis.sh`:** 1,457 lines, monolithic
- **Hard to navigate:** Many responsibilities in single script
- **Recommendation:** Split into functions or include sourced helpers

#### 9.5 Implicit Dependencies
- **Conda environment:** Script assumes `rplot` environment exists
- **Module names:** Hard-coded without documentation of why those versions
- **Recommendation:** Add `requirements.txt` and `modules.txt` for tracking

---

## 10. Recommendations by Priority

### 🔴 High Priority (Do Now)

1. **Standardize error handling:**
   - Add `set -Eeuo pipefail` to all bash scripts
   - Add trap for cleanup on error
   - Ensure consistent error propagation

2. **Add automated tests:**
   - Start with `test_prepare_combined.sh`
   - Create minimal test fixtures
   - Run in CI/CD pipeline

3. **Document module dependencies:**
   - Create `DEPENDENCIES.md` with:
     - Required HPC modules & versions
     - R packages & versions
     - System requirements

4. **Fix code duplication:**
   - Extract PIPELINE_ROOT resolution to `lib/path_resolution.sh`
   - Extract normalize functions to `lib/validation.sh`

### 🟡 Medium Priority (Next Sprint)

5. **Make QC thresholds configurable:**
   - Add `STEP1D_PCA_GENO`, `STEP1D_PCA_MIND`, `STEP1D_PCA_MAF`
   - Add LD pruning parameters to config
   - Update documentation

6. **Add sample-level QC:**
   - Create `plot_sample_metrics.R`
   - Extract per-sample stats from VCFs
   - Add `--sample-qc` mode

7. **Version checking:**
   - Verify bcftools ≥1.10, plink2 ≥2.0, R ≥4.0
   - Check R package versions
   - Warn on potential incompatibilities

8. **Improve logging:**
   - Add timing information (how long each step took)
   - Log resource usage (memory, disk)
   - Add progress indicators for long operations

### 🟢 Low Priority (Future)

9. **Parallelize per-chromosome operations:**
   - Use GNU parallel for metrics extraction
   - Parallel plot generation
   - Estimate 3-5x speedup

10. **Add resume/checkpoint system:**
    - State file tracking completed steps
    - `--resume` flag to continue from failure
    - Atomic file operations (write to `.tmp`, then `mv`)

11. **Containerization:**
    - Create Singularity definition
    - Pin all dependency versions
    - Enable portable execution

12. **Advanced features:**
    - Sample filtering by metadata
    - Customizable variant filters
    - Interactive HTML reports

---

## 11. Comparison to Coding Standards

### Pipeline-Wide Standards Compliance

| Standard | Step1D Compliance | Notes |
|----------|-------------------|-------|
| **Bash safety flags** | 🟡 Partial | Most scripts use `set -euo pipefail`, but inconsistent |
| **Logging via lib/** | ✅ Full | All scripts use `lib/logging.sh` |
| **Config via pipeline_config.sh** | ✅ Full | All params sourced from central config |
| **PIPELINE_ROOT resolution** | ✅ Full | Handles both interactive and SLURM contexts |
| **Error messages to stderr** | ✅ Full | Uses `log_error()` consistently |
| **Help text available** | ✅ Full | All scripts support `--help` |
| **Dry-run mode** | 🟡 Partial | Main script supports, utilities don't all have it |
| **Temp file cleanup** | 🟡 Partial | Done in most places, but not always on error |
| **Input validation** | 🟡 Partial | Basic checks present, could be more thorough |

### R Code Standards

| Standard | Compliance | Notes |
|----------|------------|-------|
| **Package checks** | ✅ Full | All scripts verify dependencies |
| **Error messages** | ✅ Full | Clear, actionable error messages |
| **Common utilities** | ✅ Full | Shared `common_plot_utils.R` |
| **Code style** | ✅ Full | Consistent indentation, naming |
| **Comments** | ✅ Full | Well-commented, especially complex logic |

---

## 12. Final Verdict

### Overall Assessment

**Step1D is a production-ready module with solid engineering foundations.** It demonstrates mature software development practices in most areas, with particularly strong configuration management, logging, and documentation. The recent addition of `prepare_combined_for_pca.sh` with bcftools stats generation is a significant improvement.

### Strengths Summary
- ✅ Comprehensive QC and visualization capabilities
- ✅ Flexible, configuration-driven design
- ✅ Excellent documentation (5 reference docs)
- ✅ Robust logging and error reporting
- ✅ Clean modular architecture
- ✅ Modern tooling (bcftools, PLINK2, data.table, ggplot2)

### Critical Improvements Needed
- ⚠️ Inconsistent error handling (standardize trap + `set -Eeuo pipefail`)
- ⚠️ No automated tests (blocking for production confidence)
- ⚠️ Code duplication (PIPELINE_ROOT, normalize functions)
- ⚠️ Hard-coded module names (brittleness across HPC systems)

### Score Breakdown
- **Architecture:** 4.5/5 - Clean, modular, well-integrated
- **Code Quality:** 4/5 - Good practices, but inconsistencies
- **Features:** 4.5/5 - Comprehensive, well-thought-out
- **Documentation:** 5/5 - Exceptional
- **Testing:** 1/5 - Major gap
- **Maintainability:** 4/5 - Good structure, some long scripts
- **Performance:** 4/5 - Efficient, but could parallelize better

**Final Score: 4.0/5** ⭐⭐⭐⭐

With the high-priority recommendations addressed (especially testing and error handling standardization), this module would easily achieve 4.5-5.0/5.

---

## 13. Action Items

### For the Development Team

- [ ] Create `test/step1d/` directory with unit tests
- [ ] Standardize error handling across all scripts
- [ ] Extract duplicated code to shared libraries
- [ ] Add `DEPENDENCIES.md` with version requirements
- [ ] Make QC thresholds configurable

### For Documentation

- [ ] Update `README_FIRST.md` with new `prepare_combined_for_pca.sh` integration
- [ ] Add troubleshooting section to docs
- [ ] Document expected runtimes for typical datasets
- [ ] Add migration guide for users upgrading from older versions

### For Users

- [ ] Review configuration options in `environment.template.sh`
- [ ] Test `--dry-run` mode before production runs
- [ ] Monitor log files for warnings
- [ ] Report any issues with module loading on your HPC

---

**Report Author:** AI Assistant (Cursor IDE)  
**Review Requested From:** Phu Khang Ha  
**Next Review Date:** After addressing high-priority recommendations
