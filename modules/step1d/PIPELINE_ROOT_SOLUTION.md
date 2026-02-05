# Best Solution for PIPELINE_ROOT Duplication
**Issue:** Code Duplication - PIPELINE_ROOT Resolution  
**Status:** Revised after analyzing `past_problems_and_resolutions.md`  
**Date:** 2026-02-03

---

## TL;DR - The Answer

**Use a standardized copy-paste pattern, NOT a sourced library.**

**Why?** SLURM copies scripts to `/var/spool/slurmd/`, where your library won't exist. You need PIPELINE_ROOT to find the library, but that's what you're trying to resolve (chicken-and-egg).

---

## The Problem (From Past Incidents)

From `docs/past_problems_and_resolutions.md`, PIPELINE_ROOT resolution has failed **5+ times**:

### Critical Failure: `/var/spool` Execution (§5.1, §5.5)

**What happens:**
```bash
# You run from:
/home/user/GATK_Pipeline_KH_v1/

# SLURM copies script to:
/var/spool/slurmd/job_123456/script.sh

# Your script tries:
SCRIPT_DIR="$(pwd)"  # → /var/spool/slurmd/job_123456/
PIPELINE_ROOT="${SCRIPT_DIR}/../.."  # → /var/spool/  ❌ WRONG!

source "${PIPELINE_ROOT}/lib/logging.sh"  # → File not found!
```

**Result:** Job fails with "No such file or directory"

### Why a Sourced Library Won't Work

```bash
# Idea: Create lib/path_resolution.sh and source it
source "${PIPELINE_ROOT}/lib/path_resolution.sh"  
#        ^^^^^^^^^^^^^ You don't know this yet!

# Chicken-and-egg:
# - Need PIPELINE_ROOT to find the library
# - Need the library to resolve PIPELINE_ROOT
```

---

## The Solution: Standardized Pattern v2.0

### Strategy

Instead of trying to eliminate duplication, we **standardize** it:
- ✅ Same code in every script (copy-paste from template)
- ✅ Battle-tested against all past failure modes
- ✅ Self-contained (no external dependencies)
- ✅ Version-controlled (v2.0 → know what's latest)

### The Pattern (30 lines, bulletproof)

```bash
#!/bin/bash
# PIPELINE_ROOT resolution (standardized pattern v2.0)
CONTEXT_NAME="<YOUR_SCRIPT_NAME>"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Priority 1: Explicit PIPELINE_ROOT from environment (most reliable)
if [ -n "${PIPELINE_ROOT:-}" ] && [ -d "${PIPELINE_ROOT}" ]; then
    PIPELINE_ROOT="$(cd "${PIPELINE_ROOT}" && pwd)"
    export PIPELINE_ROOT

# Priority 2: SLURM_SUBMIT_DIR (where sbatch was called from)
elif [ -n "${SLURM_SUBMIT_DIR:-}" ]; then
    if [ -f "${SLURM_SUBMIT_DIR}/config/pipeline_config.sh" ]; then
        PIPELINE_ROOT="${SLURM_SUBMIT_DIR}"
    elif [ -f "${SLURM_SUBMIT_DIR}/../config/pipeline_config.sh" ]; then
        PIPELINE_ROOT="$(cd "${SLURM_SUBMIT_DIR}/.." && pwd)"
    fi
    export PIPELINE_ROOT

# Priority 3: Script-relative discovery (interactive use)
else
    SEARCH_PATH="${SCRIPT_DIR}"
    for i in {1..3}; do
        if [ -f "${SEARCH_PATH}/config/pipeline_config.sh" ]; then
            PIPELINE_ROOT="${SEARCH_PATH}"
            break
        fi
        SEARCH_PATH="$(cd "${SEARCH_PATH}/.." && pwd)"
    done
    
    # Priority 4: $HOME fallback (last resort)
    if [ -z "${PIPELINE_ROOT:-}" ]; then
        [ -f "${HOME}/.gatk_pipeline_env" ] && source "${HOME}/.gatk_pipeline_env" 2>/dev/null || true
        PIPELINE_ROOT="${PIPELINE_ROOT:-${HOME}/GATK_Pipeline_KH_v1}"
    fi
    export PIPELINE_ROOT
fi

# Validation (fail fast with clear error)
if [ ! -f "${PIPELINE_ROOT}/config/pipeline_config.sh" ]; then
    echo "[${CONTEXT_NAME}] ❌ Invalid PIPELINE_ROOT: ${PIPELINE_ROOT}" >&2
    echo "[${CONTEXT_NAME}] Set PIPELINE_ROOT in ~/.gatk_pipeline_env or export before running" >&2
    exit 1
fi

# Resolve MODULE_DIR (avoids collision with config's SCRIPT_DIR per §5.8, §5.14)
if [[ "${SCRIPT_DIR}" =~ /modules/([^/]+) ]]; then
    MODULE_DIR="${PIPELINE_ROOT}/modules/${BASH_REMATCH[1]}"
else
    MODULE_DIR="${PIPELINE_ROOT}"
fi
```

### Search Priority Explanation

**Why this order?** Based on §5.1 in past problems:

1. **Explicit PIPELINE_ROOT** → User/admin set it (trust it)
2. **SLURM_SUBMIT_DIR** → Where job was submitted from (reliable in HPC)
3. **Script-relative** → Walk up looking for `config/pipeline_config.sh` marker
4. **$HOME fallback** → Standard location, last resort

This handles:
- ✅ Interactive execution (priority 3)
- ✅ SLURM `/var/spool` execution (priority 2)
- ✅ Moved repos (priority 1)
- ✅ Symlinked repos (priority 1)

---

## Implementation Guide

### Step 1: Create the Template

**File:** `lib/snippets/pipeline_root_resolution.snippet.sh`

(Full snippet with documentation - see IMPROVEMENT_PROPOSAL.md §1.1)

### Step 2: Update Each Script

**For:** `bin/prepare_combined_for_pca.sh`

1. **Find** the old PIPELINE_ROOT block (usually lines 1-35)
2. **Replace** with the standardized pattern
3. **Customize:**
   - Set `CONTEXT_NAME="prepare_combined_for_pca"`
   - Verify search depth is correct (usually 2-3 levels)
4. **Test:**
   ```bash
   # Interactive
   bash bin/prepare_combined_for_pca.sh /path/to/vcfs
   
   # SLURM (must export PIPELINE_ROOT!)
   sbatch --export=ALL,PIPELINE_ROOT="${PIPELINE_ROOT}" <script>
   ```

Repeat for:
- `bin/run_step1d.sh`
- `templates/master_vcf_analysis.sh`
- `templates/plink2_PCA.sh`

### Step 3: Update SLURM Wrappers

**Critical:** Always export PIPELINE_ROOT in sbatch calls:

```bash
# DO THIS:
sbatch --export=ALL,PIPELINE_ROOT="${PIPELINE_ROOT}" script.sh

# NOT THIS:
sbatch script.sh  # ❌ Script won't know PIPELINE_ROOT!
```

---

## Testing Checklist

### Test Scenarios

- [ ] **Interactive from repo root**
  ```bash
  cd /home/user/GATK_Pipeline_KH_v1
  bash modules/step1d/bin/prepare_combined_for_pca.sh --help
  ```

- [ ] **Interactive from nested directory**
  ```bash
  cd /tmp
  bash /home/user/GATK_Pipeline_KH_v1/modules/step1d/bin/prepare_combined_for_pca.sh --help
  ```

- [ ] **With PIPELINE_ROOT pre-set**
  ```bash
  export PIPELINE_ROOT=/custom/path
  bash modules/step1d/bin/prepare_combined_for_pca.sh --help
  ```

- [ ] **SLURM submission**
  ```bash
  bash modules/step1d/bin/run_step1d.sh /path/to/vcfs --PCA
  # Check logs show correct PIPELINE_ROOT
  ```

- [ ] **From `/var/spool` (simulate SLURM)**
  ```bash
  cp bin/prepare_combined_for_pca.sh /tmp/test.sh
  cd /tmp
  export SLURM_SUBMIT_DIR=/home/user/GATK_Pipeline_KH_v1
  bash test.sh --help
  # Should work via SLURM_SUBMIT_DIR
  ```

### What to Verify

1. **PIPELINE_ROOT is absolute path** (no `..` in value)
2. **Config file found** (no "missing config" errors)
3. **Libraries sourceable** (`lib/logging.sh` loads correctly)
4. **No path errors in SLURM logs** (check `/var/spool` execution)

---

## Advantages Over Library Approach

| Aspect | Sourced Library | Standardized Pattern | Winner |
|--------|----------------|---------------------|---------|
| **Works in `/var/spool`** | ❌ No (library not there) | ✅ Yes (self-contained) | Pattern |
| **No dependencies** | ❌ Needs lib/ accessible | ✅ Standalone | Pattern |
| **Handles SLURM context** | ⚠️ Complex | ✅ Built-in | Pattern |
| **Easy to audit** | ⚠️ Need to check library | ✅ Code is right there | Pattern |
| **Version tracking** | ⚠️ Library versioning | ✅ Version in comments | Pattern |
| **Code duplication** | ✅ Single source | ❌ Repeated | Library |

**Verdict:** The ~30 lines of duplication is **worth it** for reliability in HPC environments.

---

## Maintenance

### When to Update the Pattern

Update all scripts when:
- [ ] New failure mode discovered
- [ ] Search priority needs adjustment
- [ ] HPC environment changes

**How to update:**
1. Bump version number (v2.0 → v2.1)
2. Update snippet template
3. Update all scripts using the pattern
4. Add to `docs/past_problems_and_resolutions.md`

### Version History

- **v2.0** (2026-02-03) - Standardized pattern addressing all past incidents
- **v1.0** (implicit) - Original ad-hoc approaches (inconsistent)

---

## FAQ

### Q: Why not use a shell function instead of copy-paste?

**A:** Functions need to be defined before use. In `/var/spool`, you'd need to define the function in the script itself, which brings us back to copy-paste.

### Q: What if I add a 5th script that needs this?

**A:** Copy the v2.0 pattern from the snippet template. Update the version comment if you make improvements. Consider this "intentional duplication" for reliability.

### Q: Can I simplify this to just use relative paths?

**A:** No! That's what caused §5.1, §5.5 failures. Relative paths break in `/var/spool/`.

### Q: What about using `readlink -f` or `realpath`?

**A:** Good additions for Priority 1 (explicit PIPELINE_ROOT), but doesn't solve the fundamental problem that you still need to know where to look.

---

## See Also

- `docs/past_problems_and_resolutions.md` §5.1 - Pipeline root resolution failures
- `docs/past_problems_and_resolutions.md` §5.5 - SLURM /var/spool issues
- `docs/past_problems_and_resolutions.md` §5.8, §5.14 - SCRIPT_DIR collisions
- `IMPROVEMENT_PROPOSAL.md` - Full technical proposal
- `DIAGNOSTIC_REPORT.md` - Code quality analysis

---

## Quick Reference: Copy This Into Your Script

```bash
# PIPELINE_ROOT resolution (standardized pattern v2.0)
CONTEXT_NAME="<YOUR_SCRIPT_NAME>"  # ← Change this!
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [ -n "${PIPELINE_ROOT:-}" ] && [ -d "${PIPELINE_ROOT}" ]; then
    PIPELINE_ROOT="$(cd "${PIPELINE_ROOT}" && pwd)"
    export PIPELINE_ROOT
elif [ -n "${SLURM_SUBMIT_DIR:-}" ]; then
    if [ -f "${SLURM_SUBMIT_DIR}/config/pipeline_config.sh" ]; then
        PIPELINE_ROOT="${SLURM_SUBMIT_DIR}"
    elif [ -f "${SLURM_SUBMIT_DIR}/../config/pipeline_config.sh" ]; then
        PIPELINE_ROOT="$(cd "${SLURM_SUBMIT_DIR}/.." && pwd)"
    fi
    export PIPELINE_ROOT
else
    SEARCH_PATH="${SCRIPT_DIR}"
    for i in {1..3}; do
        if [ -f "${SEARCH_PATH}/config/pipeline_config.sh" ]; then
            PIPELINE_ROOT="${SEARCH_PATH}"
            break
        fi
        SEARCH_PATH="$(cd "${SEARCH_PATH}/.." && pwd)"
    done
    if [ -z "${PIPELINE_ROOT:-}" ]; then
        [ -f "${HOME}/.gatk_pipeline_env" ] && source "${HOME}/.gatk_pipeline_env" 2>/dev/null || true
        PIPELINE_ROOT="${PIPELINE_ROOT:-${HOME}/GATK_Pipeline_KH_v1}"
    fi
    export PIPELINE_ROOT
fi

if [ ! -f "${PIPELINE_ROOT}/config/pipeline_config.sh" ]; then
    echo "[${CONTEXT_NAME}] ❌ Invalid PIPELINE_ROOT: ${PIPELINE_ROOT}" >&2
    exit 1
fi
```

**Remember:** Always export PIPELINE_ROOT in sbatch!
```bash
sbatch --export=ALL,PIPELINE_ROOT="${PIPELINE_ROOT}" your_script.sh
```
