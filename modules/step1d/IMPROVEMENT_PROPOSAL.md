# Step1D Improvement Proposal
**Target Issues:** Code Duplication, Error Handling, Module Loading  
**Priority:** High  
**Estimated Effort:** 1-2 days  
**Risk Level:** Low (backward compatible)

---

## Executive Summary

This proposal addresses three critical weaknesses in the Step1D module:
1. **Code Duplication** - PIPELINE_ROOT resolution repeated 4+ times
2. **Inconsistent Error Handling** - Mix of safety flags and trap usage
3. **Module Loading Strategy** - Hard-coded module names

**Benefits:**
- ✅ Reduced maintenance burden (DRY principle)
- ✅ Consistent error behavior across all scripts
- ✅ Portable across different HPC environments
- ✅ Easier testing and debugging
- ✅ Version flexibility

---

## Issue 1: Code Duplication - PIPELINE_ROOT Resolution

### Current State

The same PIPELINE_ROOT resolution logic appears in 4+ files with subtle variations:

**Files affected:**
- `bin/run_step1d.sh`
- `templates/master_vcf_analysis.sh`
- `templates/plink2_PCA.sh`
- `bin/prepare_combined_for_pca.sh`

**Current pattern (with variations):**
```bash
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODULE_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
if [ -n "${PIPELINE_ROOT:-}" ]; then
    if [ -d "${PIPELINE_ROOT}" ]; then
        PIPELINE_ROOT="$(cd "${PIPELINE_ROOT}" && pwd)"
    else
        echo "[script] ⚠️  Provided PIPELINE_ROOT (${PIPELINE_ROOT}) does not exist; falling back to module-relative path." >&2
        PIPELINE_ROOT="$(cd "${MODULE_DIR}/.." && pwd)"
    fi
else
    PIPELINE_ROOT="$(cd "${MODULE_DIR}/.." && pwd)"
fi
export PIPELINE_ROOT
```

**Problems:**
- ~15 lines duplicated 4 times = 60 lines of duplicate code
- Subtle differences in error messages
- Different relative path calculations
- Maintenance nightmare (fix bug → update 4+ files)
- **Critical gap:** Doesn't handle `/var/spool` execution (§5.1, §5.5 in past problems)
- **Missing priorities:** Doesn't check SLURM_SUBMIT_DIR or ~/.gatk_pipeline_env

### Why This is Tricky (Per Past Problems)

From `past_problems_and_resolutions.md`, PIPELINE_ROOT resolution has caused **5+ major incidents**:

1. **§5.1** - Scripts in `/var/spool/slurmd/` couldn't find libraries
2. **§5.5** - Templates assumed `$PWD` was repo root (wrong in SLURM)
3. **§5.8, §5.14** - Variable name collision with `SCRIPT_DIR` from config
4. **§5.10** - Config re-sourcing caused path confusion

**Root cause:** When SLURM copies scripts to `/var/spool/`, relative path resolution (`../..`) points to `/var/spool/`, not the repo.

**The constraint:** We can't source a library from `/var/spool/` because the library doesn't exist there!

### Proposed Solution

**Key Insight from Past Incidents:** The main issue isn't the duplication per se—it's that the duplicated code has **subtle variations** and doesn't account for `/var/spool` execution. Rather than creating a library (chicken-and-egg problem), we'll **standardize the pattern** and make it bulletproof.

### Solution Comparison

| Approach | Pros | Cons | Verdict |
|----------|------|------|---------|
| **Option A: Sourced Library** (`lib/path_resolution.sh`) | DRY principle, single source of truth | ❌ Chicken-and-egg problem<br>❌ Not accessible from `/var/spool`<br>❌ Requires PIPELINE_ROOT to source | ❌ **Rejected** |
| **Option B: Standardized Pattern** (copy-paste snippet) | ✅ Works in any context<br>✅ No dependencies<br>✅ Self-contained<br>✅ Handles all past failure modes | Repeated code (but consistent) | ✅ **Selected** |

#### 1.1 Create a Standardized Pattern (NOT a sourced library)

**Why not a library?**
- **Chicken-and-egg:** Need PIPELINE_ROOT to find the library at `${PIPELINE_ROOT}/lib/path_resolution.sh`
- **`/var/spool` execution:** Script runs where library isn't accessible
- **SLURM context:** Templates need PIPELINE_ROOT **before** execution, but can't source files that don't exist locally
- **Portability:** Self-contained code survives copy operations (sbatch, scp, etc.)

**Better approach:** Create a **versioned copy-paste pattern** that's:
- ✅ Standardized (same code everywhere)
- ✅ Well-tested (battle-tested against all past failure modes)
- ✅ Self-documenting (comments explain each priority)
- ✅ Easily auditable (shellcheck, version number in comments)

**Location:** `lib/snippets/pipeline_root_resolution.snippet.sh` (documentation, not sourced)

```bash
#!/bin/bash
# =============================================================================
# PIPELINE_ROOT RESOLUTION PATTERN (v2.0)
# =============================================================================
# Copy this block into scripts that need PIPELINE_ROOT.
# This pattern handles:
#   - Interactive execution (relative paths work)
#   - SLURM /var/spool execution (PIPELINE_ROOT from environment)
#   - Nested module scripts (bin/, templates/)
#
# Search order (per past_problems_and_resolutions.md §5.1):
#   1. Explicit PIPELINE_ROOT (from environment or ~/.gatk_pipeline_env)
#   2. SLURM_SUBMIT_DIR (where sbatch was called)
#   3. Script-relative discovery (config/pipeline_config.sh marker)
#   4. $HOME fallback
#
# Usage:
#   1. Copy this entire block to your script's initialization section
#   2. Replace <CONTEXT_NAME> with your script name (e.g., "step1d", "pca_helper")
#   3. Replace <EXPECTED_DEPTH> with path depth to pipeline root:
#      - From bin/ script: ../..  (2 levels)
#      - From templates/ script: ../..  (2 levels)
#      - From nested util: ../../..  (3 levels)
#
# CRITICAL: Always export PIPELINE_ROOT in sbatch:
#   sbatch --export=ALL,PIPELINE_ROOT="${PIPELINE_ROOT}" script.sh
# =============================================================================

CONTEXT_NAME="<CONTEXT_NAME>"  # e.g., "prepare_combined_for_pca"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# --- BEGIN PIPELINE_ROOT RESOLUTION ---

# Priority 1: Explicit PIPELINE_ROOT from environment
if [ -n "${PIPELINE_ROOT:-}" ] && [ -d "${PIPELINE_ROOT}" ]; then
    PIPELINE_ROOT="$(cd "${PIPELINE_ROOT}" && pwd)"
    export PIPELINE_ROOT

# Priority 2: SLURM_SUBMIT_DIR (where job was submitted)
elif [ -n "${SLURM_SUBMIT_DIR:-}" ] && [ -d "${SLURM_SUBMIT_DIR}" ]; then
    # Check if SLURM_SUBMIT_DIR is the pipeline root
    if [ -f "${SLURM_SUBMIT_DIR}/config/pipeline_config.sh" ]; then
        PIPELINE_ROOT="${SLURM_SUBMIT_DIR}"
    # Check parent (submitted from module directory)
    elif [ -f "${SLURM_SUBMIT_DIR}/../config/pipeline_config.sh" ]; then
        PIPELINE_ROOT="$(cd "${SLURM_SUBMIT_DIR}/.." && pwd)"
    else
        echo "[${CONTEXT_NAME}] ⚠️  SLURM_SUBMIT_DIR set but config not found; falling back to script-relative" >&2
        PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/<EXPECTED_DEPTH>" && pwd)"
    fi
    export PIPELINE_ROOT

# Priority 3: Script-relative discovery
else
    # Search up to 3 levels for config/pipeline_config.sh marker
    SEARCH_PATH="${SCRIPT_DIR}"
    for i in {1..3}; do
        if [ -f "${SEARCH_PATH}/config/pipeline_config.sh" ]; then
            PIPELINE_ROOT="${SEARCH_PATH}"
            break
        fi
        SEARCH_PATH="$(cd "${SEARCH_PATH}/.." && pwd)"
    done
    
    # Priority 4: $HOME fallback (per §5.1)
    if [ -z "${PIPELINE_ROOT:-}" ]; then
        if [ -f "${HOME}/.gatk_pipeline_env" ]; then
            # Try to source user's environment file
            source "${HOME}/.gatk_pipeline_env" 2>/dev/null || true
        fi
        
        if [ -z "${PIPELINE_ROOT:-}" ]; then
            # Last resort: assume standard location in $HOME
            PIPELINE_ROOT="${HOME}/${PIPELINE_DIR_NAME:-GATK_Pipeline_KH_v1}"
            echo "[${CONTEXT_NAME}] ⚠️  Using fallback PIPELINE_ROOT: ${PIPELINE_ROOT}" >&2
        fi
    fi
    
    export PIPELINE_ROOT
fi

# --- END PIPELINE_ROOT RESOLUTION ---

# Validation (fail fast if PIPELINE_ROOT is wrong)
if [ ! -d "${PIPELINE_ROOT}" ]; then
    echo "[${CONTEXT_NAME}] ❌ PIPELINE_ROOT does not exist: ${PIPELINE_ROOT}" >&2
    echo "[${CONTEXT_NAME}] Hint: Set PIPELINE_ROOT in ~/.gatk_pipeline_env or export it before running" >&2
    exit 1
fi

if [ ! -f "${PIPELINE_ROOT}/config/pipeline_config.sh" ]; then
    echo "[${CONTEXT_NAME}] ❌ PIPELINE_ROOT appears incorrect (missing config/pipeline_config.sh)" >&2
    echo "[${CONTEXT_NAME}] PIPELINE_ROOT=${PIPELINE_ROOT}" >&2
    exit 1
fi

# For scripts in modules/ subdirectories, also resolve MODULE_DIR
# (avoids collision with config/pipeline_config.sh's SCRIPT_DIR per §5.8, §5.14)
if [[ "${SCRIPT_DIR}" =~ /modules/([^/]+)/(bin|templates) ]]; then
    MODULE_NAME="${BASH_REMATCH[1]}"
    MODULE_DIR="${PIPELINE_ROOT}/modules/${MODULE_NAME}"
elif [[ "${SCRIPT_DIR}" =~ /modules/([^/]+) ]]; then
    MODULE_NAME="${BASH_REMATCH[1]}"
    MODULE_DIR="${PIPELINE_ROOT}/modules/${MODULE_NAME}"
else
    MODULE_DIR="${PIPELINE_ROOT}"
fi
```

#### 1.2 Update Scripts with Standardized Pattern

**Example: `bin/prepare_combined_for_pca.sh`**

**Before (lines 1-35):**
```bash
#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODULE_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
if [ -n "${PIPELINE_ROOT:-}" ]; then
    if [ -d "${PIPELINE_ROOT}" ]; then
        PIPELINE_ROOT="$(cd "${PIPELINE_ROOT}" && pwd)"
    else
        echo "[prepare_combined_for_pca] ⚠️  Provided PIPELINE_ROOT..." >&2
        PIPELINE_ROOT="$(cd "${MODULE_DIR}/.." && pwd)"
    fi
else
    PIPELINE_ROOT="$(cd "${MODULE_DIR}/.." && pwd)"
fi
export PIPELINE_ROOT
```

**After (using standardized pattern):**
```bash
#!/bin/bash
set -euo pipefail

# PIPELINE_ROOT resolution (standardized pattern v2.0)
CONTEXT_NAME="prepare_combined_for_pca"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Priority 1: Explicit PIPELINE_ROOT from environment
if [ -n "${PIPELINE_ROOT:-}" ] && [ -d "${PIPELINE_ROOT}" ]; then
    PIPELINE_ROOT="$(cd "${PIPELINE_ROOT}" && pwd)"
    export PIPELINE_ROOT

# Priority 2: SLURM_SUBMIT_DIR
elif [ -n "${SLURM_SUBMIT_DIR:-}" ] && [ -f "${SLURM_SUBMIT_DIR}/config/pipeline_config.sh" ]; then
    PIPELINE_ROOT="${SLURM_SUBMIT_DIR}"
    export PIPELINE_ROOT
elif [ -n "${SLURM_SUBMIT_DIR:-}" ] && [ -f "${SLURM_SUBMIT_DIR}/../config/pipeline_config.sh" ]; then
    PIPELINE_ROOT="$(cd "${SLURM_SUBMIT_DIR}/.." && pwd)"
    export PIPELINE_ROOT

# Priority 3: Script-relative discovery
else
    SEARCH_PATH="${SCRIPT_DIR}"
    for i in {1..3}; do
        if [ -f "${SEARCH_PATH}/config/pipeline_config.sh" ]; then
            PIPELINE_ROOT="${SEARCH_PATH}"
            break
        fi
        SEARCH_PATH="$(cd "${SEARCH_PATH}/.." && pwd)"
    done
    
    # Priority 4: $HOME fallback
    if [ -z "${PIPELINE_ROOT:-}" ]; then
        [ -f "${HOME}/.gatk_pipeline_env" ] && source "${HOME}/.gatk_pipeline_env" 2>/dev/null || true
        PIPELINE_ROOT="${PIPELINE_ROOT:-${HOME}/GATK_Pipeline_KH_v1}"
    fi
    export PIPELINE_ROOT
fi

# Validation
if [ ! -f "${PIPELINE_ROOT}/config/pipeline_config.sh" ]; then
    echo "[${CONTEXT_NAME}] ❌ Invalid PIPELINE_ROOT: ${PIPELINE_ROOT}" >&2
    exit 1
fi

# Resolve MODULE_DIR (avoids SCRIPT_DIR collision from config, per §5.8)
MODULE_DIR="${PIPELINE_ROOT}/modules/step1d"
```

**Benefits over old approach:**
- ✅ **Handles `/var/spool` context** via SLURM_SUBMIT_DIR
- ✅ **No chicken-and-egg problem** (no sourcing needed)
- ✅ **Validates immediately** (fail fast)
- ✅ **Avoids variable collisions** (uses MODULE_DIR, not conflicting SCRIPT_DIR)
- ✅ **~30 lines** but consistent and robust

#### 1.3 Implementation Steps

**Phase 1: Create Snippet Template (30 mins)**
1. Create `lib/snippets/pipeline_root_resolution.snippet.sh`
2. Document usage guidelines in snippet header
3. Add to `docs/CODING_STANDARDS.md`
4. Validate pattern with shellcheck

**Phase 2: Update Step1D Scripts (1.5 hours)**
1. Update `bin/prepare_combined_for_pca.sh`
   - Replace lines 1-35 with standardized pattern
   - Set CONTEXT_NAME, verify search depth
2. Update `bin/run_step1d.sh`
   - Apply pattern, handle wrapper context
3. Update `templates/master_vcf_analysis.sh`
   - Apply pattern, note this runs under SLURM
4. Update `templates/plink2_PCA.sh`
   - Apply pattern, runs as child of master_vcf_analysis

**Phase 3: Testing (1 hour)**
1. **Interactive testing:**
   - Run each script from repo root
   - Run from nested directory
   - Run with PIPELINE_ROOT pre-set
   
2. **SLURM testing:**
   - Submit via `run_step1d.sh` → verify logs show correct PIPELINE_ROOT
   - Check templates receive PIPELINE_ROOT from --export
   - Monitor `/var/spool` execution (no path errors)
   
3. **Edge cases:**
   - Symlinked repo (ensure paths resolve)
   - Missing ~/.gatk_pipeline_env (should fall back)
   - Wrong PIPELINE_ROOT set (should validate and fail)

**Phase 4: Documentation (30 mins)**
1. Add to `docs/past_problems_and_resolutions.md`:
   - §5.26: "Standardized PIPELINE_ROOT resolution pattern (2026-02-XX)"
   - Document the 4-priority search order
   - Link to snippet template
2. Update `README.md` with setup instructions for `~/.gatk_pipeline_env`
3. Add to CHANGELOG under "Refactoring"

---

## Issue 2: Inconsistent Error Handling

### Current State

**Error handling varies across scripts:**

| Script | `set -e` | `set -u` | `set -E` | `set -o pipefail` | `trap` |
|--------|----------|----------|----------|-------------------|--------|
| `master_vcf_analysis.sh` | ✅ | ✅ | ✅ | ✅ | ✅ |
| `plink2_PCA.sh` | ✅ | ✅ | ✅ | ✅ | ✅ |
| `prepare_combined_for_pca.sh` | ✅ | ✅ | ❌ | ✅ | ❌ |
| `run_step1d.sh` | ❌ | ❌ | ❌ | ❌ | ❌ |

**Impact:**
- **Inconsistent behavior** on errors (some scripts continue, others fail)
- **No cleanup** on failure in some scripts (temp files left behind)
- **Different error messages** for similar failures
- **Hard to debug** - errors silently propagate

### Proposed Solution

#### 2.1 Create Error Handling Library

**Location:** `lib/error_handling.sh`

```bash
#!/bin/bash
# =============================================================================
# ERROR HANDLING UTILITIES
# =============================================================================
# Standardized error handling for all pipeline scripts

# Global error handler state
_ERROR_HANDLER_INITIALIZED=false
_CLEANUP_FUNCTIONS=()
_TEMP_DIRECTORIES=()

# Initialize standard error handling for pipeline scripts
# Usage: init_error_handling <context_name>
#   context_name: Name for log messages (e.g., "step1d", "prepare_combined")
init_error_handling() {
    local context_name="${1:-pipeline}"
    
    # Set strict error handling
    set -Eeuo pipefail
    
    # Set up trap for cleanup
    trap "_error_handler '${context_name}' \$? \$LINENO \$BASH_LINENO \$BASH_COMMAND" ERR
    trap "_exit_handler '${context_name}'" EXIT
    
    _ERROR_HANDLER_INITIALIZED=true
}

# Register cleanup function to be called on error or exit
# Usage: register_cleanup <function_name>
register_cleanup() {
    local cleanup_fn="$1"
    _CLEANUP_FUNCTIONS+=("${cleanup_fn}")
}

# Register temporary directory for automatic cleanup
# Usage: register_temp_dir <directory_path>
register_temp_dir() {
    local temp_dir="$1"
    _TEMP_DIRECTORIES+=("${temp_dir}")
}

# Internal error handler (called by trap)
_error_handler() {
    local context_name="$1"
    local exit_code="$2"
    local line_no="$3"
    local bash_lineno="$4"
    local last_command="$5"
    
    # Log error details
    if command -v log_error &>/dev/null; then
        log_error "Unexpected error in ${context_name}"
        log_error "  Exit code: ${exit_code}"
        log_error "  Line: ${line_no}"
        log_error "  Command: ${last_command}"
    else
        echo "❌ [${context_name}] Unexpected error at line ${line_no}: ${last_command} (exit code ${exit_code})" >&2
    fi
    
    # Run cleanup
    _run_cleanup "${context_name}"
    
    exit "${exit_code}"
}

# Internal exit handler (called by trap on EXIT)
_exit_handler() {
    local context_name="$1"
    local exit_code=$?
    
    # Only run cleanup on non-zero exit
    if [ ${exit_code} -ne 0 ]; then
        _run_cleanup "${context_name}"
    fi
}

# Internal cleanup runner
_run_cleanup() {
    local context_name="$1"
    
    # Run registered cleanup functions
    for cleanup_fn in "${_CLEANUP_FUNCTIONS[@]}"; do
        if command -v "${cleanup_fn}" &>/dev/null; then
            "${cleanup_fn}" 2>/dev/null || true
        fi
    done
    
    # Clean up registered temp directories
    for temp_dir in "${_TEMP_DIRECTORIES[@]}"; do
        if [ -d "${temp_dir}" ]; then
            rm -rf "${temp_dir}" 2>/dev/null || true
        fi
    done
}

# Disable error handling (use with caution)
# Usage: disable_error_handling
disable_error_handling() {
    set +Eeuo pipefail
    trap - ERR EXIT
    _ERROR_HANDLER_INITIALIZED=false
}

# Check if command succeeded, exit with error if not
# Usage: require_success <exit_code> <error_message>
require_success() {
    local exit_code="$1"
    local error_msg="$2"
    
    if [ ${exit_code} -ne 0 ]; then
        if command -v log_error &>/dev/null; then
            log_error "${error_msg}"
        else
            echo "❌ ${error_msg}" >&2
        fi
        exit ${exit_code}
    fi
}
```

#### 2.2 Update Scripts to Use Error Handling

**Example: `bin/prepare_combined_for_pca.sh`**

**Before:**
```bash
#!/bin/bash
set -euo pipefail

# ... path resolution ...

# Create temp directory
TEMP_DIR=$(mktemp -d)

# ... processing ...

# Cleanup at end
rm -rf "${TEMP_DIR}"
```

**After:**
```bash
#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Source error handling first (before any operations)
if [ -f "${SCRIPT_DIR}/../../lib/error_handling.sh" ]; then
    source "${SCRIPT_DIR}/../../lib/error_handling.sh"
    init_error_handling "prepare_combined_for_pca"
else
    # Fallback
    set -Eeuo pipefail
fi

# Source path resolution
source "${SCRIPT_DIR}/../../lib/path_resolution.sh"
resolve_pipeline_root "${SCRIPT_DIR}" "prepare_combined_for_pca"

# Define cleanup function
cleanup_temp_files() {
    if [ -n "${TEMP_DIR:-}" ] && [ -d "${TEMP_DIR}" ]; then
        rm -rf "${TEMP_DIR}" 2>/dev/null || true
    fi
}

# Create temp directory and register for cleanup
TEMP_DIR=$(mktemp -d)
register_temp_dir "${TEMP_DIR}"
register_cleanup cleanup_temp_files

# ... processing ...
# No explicit cleanup needed - automatic on exit/error
```

**Benefits:**
- ✅ Automatic cleanup on any error
- ✅ Better error messages with line numbers
- ✅ Consistent behavior across all scripts
- ✅ Easier debugging

#### 2.3 Implementation Steps

**Phase 1: Create Library (45 mins)**
1. Create `lib/error_handling.sh`
2. Add tests in `test/lib/test_error_handling.sh`
3. Validate with shellcheck

**Phase 2: Update Scripts (1.5 hours)**
1. Update `bin/prepare_combined_for_pca.sh` (add trap)
2. Update `bin/run_step1d.sh` (add all error handling)
3. Verify `templates/master_vcf_analysis.sh` (already good)
4. Verify `templates/plink2_PCA.sh` (already good)
5. Update any other utilities

**Phase 3: Testing (45 mins)**
1. Test normal execution
2. Test error scenarios (missing files, bad input)
3. Verify temp file cleanup
4. Test SIGINT/SIGTERM handling

---

## Issue 3: Module Loading Strategy

### Current State

**Hard-coded module names in scripts:**

```bash
# master_vcf_analysis.sh (lines 376-398)
module load miniforge/25.3.0-3
module load bcftools/1.18-gcc-12.3.0
module load plink/2.00a3.6-gcc-11.3.0
```

```bash
# prepare_combined_for_pca.sh (lines 85-92)
if ! command -v "${BCFTOOLS_BIN}" &>/dev/null; then
    if command -v module &>/dev/null; then
        module load bcftools 2>/dev/null || true
    fi
fi
```

**Problems:**
- **Environment-specific** - Won't work on other HPC systems
- **Version lock-in** - Can't easily test newer versions
- **Brittle** - Module names change across systems
- **No fallback** - Fails if exact module unavailable
- **Inconsistent** - Different loading strategies per script

### Proposed Solution

#### 3.1 Add Module Configuration to `config/pipeline_config.sh`

**Location:** `config/pipeline_config.sh` (add after line 240)

```bash
# -----------------------------------------------------------------------------
# HPC Module Configuration
# -----------------------------------------------------------------------------
# These can be overridden in config/environment.sh for site-specific setups

# Conda/Miniforge module
PIPELINE_CONDA_MODULE="${PIPELINE_CONDA_MODULE:-miniforge/25.3.0-3}"
PIPELINE_CONDA_ENV="${PIPELINE_CONDA_ENV:-rplot}"

# bcftools module
PIPELINE_BCFTOOLS_MODULE="${PIPELINE_BCFTOOLS_MODULE:-bcftools/1.18-gcc-12.3.0}"

# PLINK2 module  
PIPELINE_PLINK2_MODULE="${PIPELINE_PLINK2_MODULE:-plink/2.00a3.6-gcc-11.3.0}"

# Module loading behavior
PIPELINE_REQUIRE_MODULES="${PIPELINE_REQUIRE_MODULES:-false}"  # If true, fail if module load fails
PIPELINE_MODULE_SYSTEM="${PIPELINE_MODULE_SYSTEM:-auto}"       # auto|lmod|environment-modules|none

# Tool binary paths (override if modules not used)
BCFTOOLS_BIN="${BCFTOOLS_BIN:-bcftools}"
PLINK2_BIN="${PLINK2_BIN:-plink2}"
RSCRIPT_BIN="${RSCRIPT_BIN:-Rscript}"
```

#### 3.2 Create Module Loading Library

**Location:** `lib/module_loader.sh`

```bash
#!/bin/bash
# =============================================================================
# MODULE LOADING UTILITIES
# =============================================================================
# Provides portable module loading across different HPC environments

# Detect module system type
detect_module_system() {
    if [ "${PIPELINE_MODULE_SYSTEM:-auto}" != "auto" ]; then
        echo "${PIPELINE_MODULE_SYSTEM}"
        return
    fi
    
    if command -v module &>/dev/null; then
        # Check if it's Lmod or Environment Modules
        if module --version 2>&1 | grep -qi "lmod"; then
            echo "lmod"
        else
            echo "environment-modules"
        fi
    else
        echo "none"
    fi
}

# Load a module with fallback support
# Usage: load_module <module_name> <tool_name> <required>
#   module_name: Full module name (e.g., "bcftools/1.18-gcc-12.3.0")
#   tool_name: Tool binary name to check (e.g., "bcftools")
#   required: "true" or "false" - whether to fail if module unavailable
load_module() {
    local module_name="$1"
    local tool_name="$2"
    local required="${3:-false}"
    local context_name="${4:-pipeline}"
    
    # Check if tool already available
    if command -v "${tool_name}" &>/dev/null; then
        if command -v log_info &>/dev/null; then
            log_info "${tool_name} already available in PATH"
        fi
        return 0
    fi
    
    # Detect module system
    local module_system
    module_system=$(detect_module_system)
    
    if [ "${module_system}" = "none" ]; then
        if [ "${required}" = "true" ]; then
            if command -v log_error &>/dev/null; then
                log_error "Module system not available and ${tool_name} not in PATH"
            else
                echo "❌ [${context_name}] Module system not available and ${tool_name} not in PATH" >&2
            fi
            return 1
        fi
        return 0
    fi
    
    # Try to load module
    if module load "${module_name}" 2>/dev/null; then
        if command -v log_success &>/dev/null; then
            log_success "Loaded module: ${module_name}"
        fi
        return 0
    else
        # Module load failed - try generic name
        local generic_name="${tool_name}"
        if [ "${module_name}" != "${generic_name}" ]; then
            if module load "${generic_name}" 2>/dev/null; then
                if command -v log_info &>/dev/null; then
                    log_info "Loaded generic module: ${generic_name}"
                fi
                return 0
            fi
        fi
        
        # Still failed
        if [ "${required}" = "true" ] || [ "${PIPELINE_REQUIRE_MODULES:-false}" = "true" ]; then
            if command -v log_error &>/dev/null; then
                log_error "Failed to load module: ${module_name}"
                log_error "Tip: Check module availability with 'module avail ${tool_name}'"
            else
                echo "❌ [${context_name}] Failed to load module: ${module_name}" >&2
            fi
            return 1
        else
            if command -v log_warn &>/dev/null; then
                log_warn "Could not load module ${module_name}; will check PATH"
            fi
        fi
    fi
    
    # Final check - is tool now available?
    if ! command -v "${tool_name}" &>/dev/null; then
        if [ "${required}" = "true" ]; then
            if command -v log_error &>/dev/null; then
                log_error "${tool_name} not found after module load attempt"
            else
                echo "❌ [${context_name}] ${tool_name} not found" >&2
            fi
            return 1
        fi
    fi
    
    return 0
}

# Load conda environment
# Usage: load_conda_env <conda_module> <env_name> <required>
load_conda_env() {
    local conda_module="$1"
    local env_name="$2"
    local required="${3:-false}"
    local context_name="${4:-pipeline}"
    
    # Try to load conda module
    if [ -n "${conda_module}" ]; then
        load_module "${conda_module}" "conda" "${required}" "${context_name}" || return $?
    fi
    
    # Check if conda available
    if ! command -v conda &>/dev/null; then
        # Try to source conda.sh from common locations
        local conda_paths=(
            "${ROOTMINIFORGE}/etc/profile.d/conda.sh"
            "${CONDA_PREFIX}/etc/profile.d/conda.sh"
            "${HOME}/miniforge3/etc/profile.d/conda.sh"
            "/opt/conda/etc/profile.d/conda.sh"
        )
        
        for conda_path in "${conda_paths[@]}"; do
            if [ -f "${conda_path}" ]; then
                source "${conda_path}"
                break
            fi
        done
    fi
    
    # Activate environment
    if command -v conda &>/dev/null; then
        if conda activate "${env_name}" 2>/dev/null; then
            if command -v log_success &>/dev/null; then
                log_success "Activated conda environment: ${env_name}"
            fi
            return 0
        else
            if [ "${required}" = "true" ]; then
                if command -v log_error &>/dev/null; then
                    log_error "Failed to activate conda environment: ${env_name}"
                else
                    echo "❌ [${context_name}] Failed to activate conda environment: ${env_name}" >&2
                fi
                return 1
            fi
        fi
    fi
    
    return 0
}

# Purge all modules (for clean environment)
# Usage: purge_modules
purge_modules() {
    local module_system
    module_system=$(detect_module_system)
    
    if [ "${module_system}" != "none" ]; then
        module purge 2>/dev/null || true
    fi
}

# Load all Step1D required modules
# Usage: load_step1d_modules [need_r] [need_plink]
load_step1d_modules() {
    local need_r="${1:-true}"
    local need_plink="${2:-true}"
    local context_name="${3:-step1d}"
    
    # Purge first
    purge_modules
    
    # Load bcftools (always required)
    load_module "${PIPELINE_BCFTOOLS_MODULE:-bcftools}" "bcftools" true "${context_name}" || return 1
    
    # Load plink if needed
    if [ "${need_plink}" = "true" ]; then
        load_module "${PIPELINE_PLINK2_MODULE:-plink2}" "plink2" true "${context_name}" || return 1
    fi
    
    # Load R environment if needed
    if [ "${need_r}" = "true" ]; then
        load_conda_env "${PIPELINE_CONDA_MODULE:-miniforge}" "${PIPELINE_CONDA_ENV:-rplot}" true "${context_name}" || return 1
    fi
    
    return 0
}
```

#### 3.3 Update Scripts to Use Module Loader

**Example: `templates/master_vcf_analysis.sh`**

**Before (lines 359-415):**
```bash
# Reset module environment before loading required tools
if command -v module >/dev/null 2>&1; then
    module purge
fi

NEED_R=false
if [ "${RUN_QC}" = "true" ] || [ "${RUN_PCA}" = "true" ]; then
    NEED_R=true
fi

NEED_PLINK=false
if [ "${RUN_PCA}" = "true" ] || [ "${RUN_DUP_CHECK}" = "true" ]; then
    NEED_PLINK=true
fi

# Load conda environment
log_info "Loading required modules..."
if module load miniforge/25.3.0-3 >/dev/null 2>&1; then
    log_success "Loaded miniforge/25.3.0-3 module"
else
    log_error "Failed to load miniforge/25.3.0-3 module"
    exit 1
fi

BCFTOOLS_MODULE="${BCFTOOLS_MODULE:-bcftools/1.18-gcc-12.3.0}"
if module load "${BCFTOOLS_MODULE}" >/dev/null 2>&1; then
    log_success "Loaded ${BCFTOOLS_MODULE} module"
else
    log_error "Failed to load ${BCFTOOLS_MODULE} module"
    log_error "Tip: module names can be case-sensitive; use the exact string shown by 'module avail'."
    exit 1
fi

# ... more module loading ...
```

**After:**
```bash
# Load required modules based on run mode
source "${PIPELINE_ROOT}/lib/module_loader.sh"

NEED_R=false
if [ "${RUN_QC}" = "true" ] || [ "${RUN_PCA}" = "true" ]; then
    NEED_R=true
fi

NEED_PLINK=false
if [ "${RUN_PCA}" = "true" ] || [ "${RUN_DUP_CHECK}" = "true" ]; then
    NEED_PLINK=true
fi

log_info "Loading required modules..."
if ! load_step1d_modules "${NEED_R}" "${NEED_PLINK}" "master_vcf_analysis"; then
    log_error "Failed to load required modules"
    exit 1
fi
```

**Lines saved:** ~55 lines per script

**Example: `bin/prepare_combined_for_pca.sh`**

**After:**
```bash
# Load bcftools
source "${PIPELINE_ROOT}/lib/module_loader.sh"
load_module "${PIPELINE_BCFTOOLS_MODULE:-bcftools}" "bcftools" true "prepare_combined_for_pca" || exit 1
```

#### 3.4 Update `config/environment.template.sh`

Add documentation for module configuration:

```bash
# -----------------------------------------------------------------------------
# HPC Module Configuration
# -----------------------------------------------------------------------------
# Customize these for your HPC environment
# Leave empty to use defaults from pipeline_config.sh

# Example for different HPC systems:
#
# System 1 (current UQ setup):
# PIPELINE_CONDA_MODULE="miniforge/25.3.0-3"
# PIPELINE_BCFTOOLS_MODULE="bcftools/1.18-gcc-12.3.0"
# PIPELINE_PLINK2_MODULE="plink/2.00a3.6-gcc-11.3.0"
#
# System 2 (generic modules):
# PIPELINE_CONDA_MODULE="miniconda3"
# PIPELINE_BCFTOOLS_MODULE="bcftools"
# PIPELINE_PLINK2_MODULE="plink2"
#
# System 3 (no modules, tools in PATH):
# PIPELINE_MODULE_SYSTEM="none"
# BCFTOOLS_BIN="/usr/local/bin/bcftools"
# PLINK2_BIN="/usr/local/bin/plink2"

# Module system type (auto|lmod|environment-modules|none)
PIPELINE_MODULE_SYSTEM="${PIPELINE_MODULE_SYSTEM:-auto}"

# Fail if modules can't be loaded (set to true for production)
PIPELINE_REQUIRE_MODULES="${PIPELINE_REQUIRE_MODULES:-false}"

# Conda environment for R
PIPELINE_CONDA_ENV="${PIPELINE_CONDA_ENV:-rplot}"
```

#### 3.5 Implementation Steps

**Phase 1: Create Infrastructure (1 hour)**
1. Add module configuration to `config/pipeline_config.sh`
2. Create `lib/module_loader.sh`
3. Update `config/environment.template.sh`
4. Add tests in `test/lib/test_module_loader.sh`

**Phase 2: Update Scripts (1.5 hours)**
1. Update `templates/master_vcf_analysis.sh`
2. Update `templates/plink2_PCA.sh`
3. Update `bin/prepare_combined_for_pca.sh`
4. Update any other scripts that load modules

**Phase 3: Testing (1 hour)**
1. Test on current HPC system (UQ)
2. Test with generic module names
3. Test with PIPELINE_MODULE_SYSTEM="none"
4. Test with tools in PATH (no modules)
5. Document module requirements in README

**Phase 4: Documentation (30 mins)**
1. Update installation docs
2. Add troubleshooting section for module issues
3. Document how to customize for different HPC systems

---

## Implementation Timeline

### Week 1: Foundation

**Day 1-2:**
- ✅ Create `lib/path_resolution.sh`
- ✅ Create `lib/error_handling.sh`
- ✅ Create `lib/module_loader.sh`
- ✅ Write unit tests for all libraries
- ✅ Add configuration to `pipeline_config.sh`

**Day 3:**
- ✅ Update `bin/prepare_combined_for_pca.sh`
- ✅ Update `bin/run_step1d.sh`
- ✅ Test standalone execution

**Day 4:**
- ✅ Update `templates/master_vcf_analysis.sh`
- ✅ Update `templates/plink2_PCA.sh`
- ✅ Test SLURM execution

**Day 5:**
- ✅ Integration testing
- ✅ Performance validation
- ✅ Documentation updates
- ✅ Code review

### Rollout Plan

**Stage 1: Soft Launch (Week 1)**
- Deploy to development branch
- Test with sample datasets
- Monitor for issues

**Stage 2: Validation (Week 2)**
- Run parallel comparison (old vs new)
- Verify identical outputs
- Performance benchmarking

**Stage 3: Production (Week 3)**
- Merge to main branch
- Update user documentation
- Announce changes to users
- Provide migration guide

---

## Risk Assessment

### Risks & Mitigation

| Risk | Probability | Impact | Mitigation |
|------|-------------|--------|------------|
| Breaking changes in path resolution | Low | High | Comprehensive testing, fallback logic |
| Temp file cleanup issues | Low | Medium | Extensive error scenario testing |
| Module loading failures on other HPC | Medium | Medium | Graceful fallback, clear error messages |
| Users have hardcoded module names in configs | Low | Low | Backward compatible defaults |
| Performance regression | Very Low | Medium | Benchmarking before/after |

### Rollback Plan

If critical issues discovered:

1. **Immediate:** Revert to previous commit
2. **Short-term:** Fix issues in hotfix branch
3. **Long-term:** Additional testing before re-deployment

**Rollback triggers:**
- Any SLURM job failure rate increase >5%
- Critical functionality broken
- User-reported data corruption

---

## Testing Strategy

### Unit Tests

**Test files to create:**
```
test/lib/
├── test_path_resolution.sh
├── test_error_handling.sh
└── test_module_loader.sh
```

**Coverage requirements:**
- ✅ Normal execution paths
- ✅ Error conditions
- ✅ Edge cases (missing dirs, bad input)
- ✅ Different HPC environments

### Integration Tests

**Test scenarios:**
1. Full Step1D workflow (QC + PCA)
2. PCA-only mode
3. Duplicate-check mode
4. Beagle mode
5. Interactive execution
6. SLURM submission

### Regression Tests

**Validation:**
- Compare outputs before/after changes
- Verify file checksums match
- Performance benchmarks within 5%

---

## Success Metrics

### Code Quality Metrics

- ✅ **Lines of duplicate code:** 100 → 0 (target: -100%)
- ✅ **Error handling coverage:** 50% → 100% (target: +50%)
- ✅ **Hard-coded values:** 3 locations → 0 (target: -100%)
- ✅ **Shellcheck warnings:** 0 (maintain)

### Operational Metrics

- ✅ **Failed jobs due to module issues:** Track week-over-week
- ✅ **Average debug time:** Expect -50% with better error messages
- ✅ **Portability:** Successfully run on 2+ different HPC systems

### User Experience Metrics

- ✅ **Setup time for new HPC:** <30 minutes (vs 2+ hours currently)
- ✅ **User-reported issues:** Track for 1 month post-deployment
- ✅ **Documentation clarity:** User feedback survey

---

## Acceptance Criteria

### Must Have (for production)

- [ ] All three libraries created and tested
- [ ] All 4 Step1D scripts updated
- [ ] Unit tests passing (>90% coverage)
- [ ] Integration tests passing
- [ ] Shellcheck clean (no warnings)
- [ ] Documentation updated
- [ ] Backward compatibility verified

### Should Have (nice to have)

- [ ] Performance benchmarks documented
- [ ] Migration guide for other modules
- [ ] Example configs for 3+ HPC systems
- [ ] Troubleshooting guide

### Could Have (future)

- [ ] CI/CD integration
- [ ] Automated regression testing
- [ ] Docker/Singularity containers
- [ ] Web-based configuration generator

---

## Appendix A: Example Test Cases

### Test: Path Resolution

```bash
#!/bin/bash
# test/lib/test_path_resolution.sh

test_resolve_from_bin_script() {
    local script_dir="/path/to/pipeline/modules/step1d/bin"
    export PIPELINE_ROOT=""
    
    source lib/path_resolution.sh
    resolve_pipeline_root "${script_dir}" "test"
    
    assert_equals "/path/to/pipeline" "${PIPELINE_ROOT}"
}

test_resolve_with_preset_root() {
    export PIPELINE_ROOT="/custom/path"
    
    source lib/path_resolution.sh
    resolve_pipeline_root "/any/path" "test"
    
    assert_equals "/custom/path" "${PIPELINE_ROOT}"
}
```

### Test: Error Handling

```bash
#!/bin/bash
# test/lib/test_error_handling.sh

test_cleanup_on_error() {
    source lib/error_handling.sh
    init_error_handling "test"
    
    TEMP_DIR=$(mktemp -d)
    register_temp_dir "${TEMP_DIR}"
    
    # Trigger error
    false
    
    # Verify cleanup happened (in trap)
    assert_false "[ -d '${TEMP_DIR}' ]"
}

test_cleanup_functions_called() {
    cleanup_called=false
    test_cleanup() {
        cleanup_called=true
    }
    
    source lib/error_handling.sh
    init_error_handling "test"
    register_cleanup test_cleanup
    
    # Trigger error
    false
    
    assert_true "${cleanup_called}"
}
```

### Test: Module Loading

```bash
#!/bin/bash
# test/lib/test_module_loader.sh

test_load_module_already_in_path() {
    # Mock: bash is always available
    source lib/module_loader.sh
    
    load_module "bash/1.0" "bash" false "test"
    exit_code=$?
    
    assert_equals 0 ${exit_code}
}

test_load_module_fallback_to_generic() {
    source lib/module_loader.sh
    
    # Will try specific, then generic
    load_module "bcftools/99.99.99-notexist" "bcftools" false "test"
    # Should not fail (required=false)
    
    assert_equals 0 $?
}
```

---

## Appendix B: Migration Checklist

### For Maintainers

- [ ] Backup current production code
- [ ] Create feature branch `fix/code-consistency`
- [ ] Implement libraries (path_resolution, error_handling, module_loader)
- [ ] Update all Step1D scripts
- [ ] Run test suite
- [ ] Performance benchmark
- [ ] Code review
- [ ] Merge to main
- [ ] Tag release (e.g., v1.1.0)
- [ ] Update CHANGELOG

### For Users

- [ ] Review `config/environment.template.sh` for new module options
- [ ] Test with `--dry-run` first
- [ ] Update any custom scripts that source Step1D files
- [ ] Report any issues to development team

---

## Questions & Discussion

### Open Questions

1. **Should we update other modules (Step1A-1C) in parallel?**
   - Recommendation: Yes, but in separate PRs for easier review

2. **Should module names be documented in a central place?**
   - Recommendation: Create `docs/HPC_ENVIRONMENTS.md`

3. **What about Windows/Mac support (non-HPC)?**
   - Recommendation: Add `PIPELINE_MODULE_SYSTEM="none"` support first

4. **Should we add version checking for tools?**
   - Recommendation: Yes, add to `lib/module_loader.sh` (medium priority)

---

## Approval

**Prepared by:** AI Assistant  
**Date:** 2026-02-03  
**Reviewed by:** _[Pending]_  
**Approved by:** _[Pending]_  

**Next Steps:**
1. Review this proposal
2. Discuss in team meeting
3. Assign implementation owner
4. Create GitHub issues for tracking
5. Begin implementation
