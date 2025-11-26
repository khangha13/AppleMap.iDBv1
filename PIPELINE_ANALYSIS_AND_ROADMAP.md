# GATK Pipeline KH v1 - AI Agent Implementation Guide

**Purpose**: This document provides clear instructions for AI agents working on improving the GATK Pipeline KH v1 codebase.

**Author**: Phu Khang HA  
**Last Updated**: 2025-01-11

---

## 🎯 OBJECTIVE

Transform the GATK Pipeline KH v1 into a consistent, easy-to-use pipeline with:
- Single entry point supporting all 4 steps (1A, 1B, 1C, 1D)
- Both CLI and interactive interfaces for each step
- Clean directory structure with no redundant files
- Complete documentation

---

## 📊 CURRENT STATE ANALYSIS

### **What Works Well**
- Modular architecture: Each step (1A-1D) has consistent structure (`bin/`, `lib/`, `templates/`)
- Professional logging system (`lib/logging.sh`) with SLURM integration
- Centralized configuration (`config/pipeline_config.sh`)
- SLURM job management with dynamic script generation

### **What Needs Fixing**

1. **Main entry point incomplete**
   - `GATK_pipeline.sh` only handles Steps 1A and 1B
   - Steps 1C (Beagle) and 1D (VCF Analysis) exist but are not integrated
   - No CLI flags, only interactive prompts

2. **Missing interactive wrappers**
   - Only `wrappers/interactive/step1d_interactive.sh` exists
   - Need: `step1a_interactive.sh`, `step1b_interactive.sh`, `step1c_interactive.sh`

3. **Redundant files**
   - `Original_scripts/` contains old monolithic scripts (replaced by modules)
   - `redundant/` folder with superseded configs
   - `.DS_Store` system files

4. **Documentation gaps**
   - README.md only mentions Steps 1A and 1B
   - Steps 1C and 1D not documented in main README

---

## 📁 TARGET DIRECTORY STRUCTURE

**IMPORTANT**: This is the target structure. All work should move toward this layout.

```
GATK_Pipeline_KH_v1/
│
├── bin/                       # All user-facing executables
│   └── gatk_pipeline.sh       # single entry point (CLI)
│
├── modules/                   # step-specific logic (unchanged)
│   ├── step1a/
│   ├── step1b/
│   ├── step1c/
│   └── step1d/
│
├── wrappers/
│   ├── interactive/           # per-step interactive helpers
│   │   ├── step1a_interactive.sh
│   │   ├── step1b_interactive.sh
│   │   ├── step1c_interactive.sh
│   │   └── step1d_interactive.sh
│   └── sbatch/                # minimal sbatch submitters
│       ├── step1a_submit.sh
│       ├── step1b_submit.sh
│       ├── step1c_submit.sh
│       └── step1d_submit.sh
│
├── config/
│   ├── env.template.sh        # user fills then sources
│   └── pipeline_config.sh
│
├── lib/                       # common bash libraries
│   ├── logging.sh
│   ├── validation.sh
│   ├── slurm.sh
│   └── pipeline_common.sh
│
├── docs/                      # markdown docs only
├── test/                      # automated test data + bats tests
└── archive/                   # moved original & redundant scripts
```

### **Key Structural Changes**

1. **`GATK_pipeline.sh` → `bin/gatk_pipeline.sh`**
   - Move main entry point to `bin/` directory
   - Make executable (chmod +x)
   - Standard location for user-facing scripts

2. **`archive/` directory**
   - New location for redundant/original scripts
   - Preserves history without cluttering main directory

---

## 🗂️ FILES TO ARCHIVE

**ACTION**: Move these to `archive/` directory before starting improvements.

### **Safe to Archive Immediately**

1. **`Original_scripts/`** (entire directory)
   - `Apple_GATK_Pipeline_1a_CallVariantsPerSampleKHv2.sh` → Replaced by `modules/step1a/`
   - `Apple_GATK_Pipeline_1B_CombineGVCFs.sh` → Replaced by `modules/step1b/`
   - `Wrapper_GATK.sh` → Replaced by `GATK_pipeline.sh`
   - `Conda wrapper script` → Functionality integrated into modules

2. **`redundant/`** (entire directory)
   - `step1a_config.sh` → Replaced by `config/pipeline_config.sh`
   - `step1b_config.sh` → Replaced by `config/pipeline_config.sh`
   - `step1d_samtools_pull_bedfile.sh` → Functionality moved to modules

3. **System Files**
   - All `.DS_Store` files → Delete and add to `.gitignore`

### **Review Before Archiving**

1. **`test/` directory** - Keep for now, but consider migrating to `bats-core` framework
2. **`utils/` directory** - Check if R scripts duplicate `modules/step1d/Rscripts/` and consolidate if needed

---

## 🔧 IMPLEMENTATION TASKS

### **TASK 0: Repository Cleanup**

**Goal**: Clean up redundant files and create proper directory structure.

**Steps**:
1. Delete all `.DS_Store` files:
   ```bash
   find . -name .DS_Store -delete
   ```

2. Add `.DS_Store` to `.gitignore`:
   ```bash
   echo ".DS_Store" >> .gitignore
   ```

3. Create `archive/` directory:
   ```bash
   mkdir -p archive
   ```

4. Move redundant directories:
   ```bash
   mv Original_scripts archive/original_scripts
   mv redundant archive/redundant
   ```

5. Create `bin/` directory:
   ```bash
   mkdir -p bin
   ```

6. Add code quality configuration files:
   - Create `.editorconfig` with standard bash formatting rules
   - Create `.shellcheckrc` with: `shellcheck --shell=bash --exclude=SC1090,SC1091`

**Verification**: Check that `archive/`, `bin/` exist, `.DS_Store` is in `.gitignore`, and redundant files are moved.

---

### **TASK 1: Create Unified Entry Point**

**Goal**: Transform `GATK_pipeline.sh` into `bin/gatk_pipeline.sh` with CLI support and full step integration.

**Steps**:

1. **Read and understand current `GATK_pipeline.sh`**
   - Note how it currently handles Steps 1A and 1B
   - Understand the interactive prompt flow
   - Identify functions that need to be preserved

2. **Add CLI argument parsing**
   - Use `getopt` or manual parsing (prefer `getopt` if available)
   - Support these flags:
     - `-d, --dataset DATASET` (required when not in interactive mode)
     - `-s, --step {1a|1b|1c|1d|full}` (default: auto-detect)
     - `-i, --interactive` (force interactive mode)
     - `--dry-run` (preview without execution)
     - `-h, --help` (show usage and exit)
     - `-q, --quiet` (suppress banners and non-essential output)

3. **Integrate Steps 1C and 1D**
   - Add functions `run_step1c()` and `run_step1d()` to `lib/pipeline_common.sh`
   - These should call `wrappers/sbatch/step1c_submit.sh` and `step1d_submit.sh` respectively
   - Follow the same pattern as `run_step1a()` and `run_step1b()`

4. **Add status checking for Steps 1C and 1D**
   - In `lib/validation.sh`, add:
     - `check_step1c_status(rdm_base_path)` - Check for Beagle output files
     - `check_step1d_status(rdm_base_path)` - Check for VCF analysis outputs
   - Update `analyze_pipeline_status()` in `lib/pipeline_common.sh` to include all 4 steps

5. **Update main workflow logic**
   - Modify `route_pipeline_execution()` to handle Steps 1C and 1D
   - Add handlers: `handle_step1c_needed()`, `handle_step1d_needed()`
   - Update `run_complete_pipeline()` to optionally include Steps 1C and 1D

6. **Maintain backward compatibility**
   - If no CLI arguments provided, default to interactive mode
   - Preserve all existing interactive prompts
   - Ensure existing workflows continue to work

7. **Move and rename**
   - Copy `GATK_pipeline.sh` to `bin/gatk_pipeline.sh`
   - Make executable: `chmod +x bin/gatk_pipeline.sh`
   - Update all internal path references (SCRIPT_DIR calculations)
   - Keep original `GATK_pipeline.sh` temporarily for reference, then remove after verification

**Verification**:
- `bin/gatk_pipeline.sh --help` shows usage
- `bin/gatk_pipeline.sh -d test_dataset -s 1a` runs Step 1A non-interactively
- `bin/gatk_pipeline.sh` (no args) enters interactive mode
- All 4 steps can be executed via CLI flags

---

### **TASK 2: Create Missing Interactive Wrappers**

**Goal**: Create interactive wrapper scripts for Steps 1A, 1B, and 1C following the pattern of `step1d_interactive.sh`.

**Reference**: Study `wrappers/interactive/step1d_interactive.sh` to understand the pattern.

**Steps for each wrapper**:

#### **2.1 Create `wrappers/interactive/step1a_interactive.sh`**

**Requirements**:
- Prompt user for dataset name
- Prompt user for RDM base path (or infer from dataset name)
- Validate path exists
- Show sample list from `1.FASTQ/` directory
- Show completion status (which samples are done)
- Confirm before calling `wrappers/sbatch/step1a_submit.sh`
- Use consistent formatting and error messages

**Implementation pattern**:
```bash
#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

# Source common functions
source "${PIPELINE_ROOT}/lib/logging.sh"
source "${PIPELINE_ROOT}/lib/validation.sh"
source "${PIPELINE_ROOT}/config/pipeline_config.sh"

# Prompt functions
# Validation functions
# Main execution logic
```

#### **2.2 Create `wrappers/interactive/step1b_interactive.sh`**

**Requirements**:
- Prompt for dataset name and RDM base path
- Show chromosome list (from reference genome)
- Show completion status (which chromosomes are done)
- Validate Step 1A outputs exist (`5.Individual_VCF/`)
- Confirm before calling `wrappers/sbatch/step1b_submit.sh`

#### **2.3 Create `wrappers/interactive/step1c_interactive.sh`**

**Requirements**:
- Prompt for dataset name and RDM base path
- Prompt for gene map file (optional, with validation)
- Show VCF manifest (list of VCFs to process)
- Show output directory
- Confirm before calling `wrappers/sbatch/step1c_submit.sh`

#### **2.4 Standardize All Interactive Scripts**

**Create common functions** (consider adding to `lib/pipeline_common.sh` or a new `lib/interactive_common.sh`):
- `prompt_for_dataset()` - Standardized dataset name prompt
- `prompt_for_path()` - Standardized path prompt with validation
- `confirm_action(message)` - Standardized yes/no confirmation
- `show_step_summary()` - Standardized summary display

**Apply to all interactive scripts**:
- Use common functions where possible
- Consistent error messages
- Consistent formatting (use same emoji/icons)
- Unified help/usage output

**Verification**:
- All 4 interactive scripts exist in `wrappers/interactive/`
- Each script can be run standalone
- Each script validates inputs before proceeding
- Consistent user experience across all scripts

---

### **TASK 3: Update Status Checking Functions**

**Goal**: Add status checking for Steps 1C and 1D to enable automatic workflow detection.

**Steps**:

1. **Add `check_step1c_status()` to `lib/validation.sh`**
   - Check for Beagle output files in expected location
   - Return: "Complete", "Partial", "Not Started"
   - Follow pattern of `check_step1a_status()` and `check_step1b_status()`
   - Check for expected file patterns (consult `modules/step1c/` to understand outputs)

2. **Add `check_step1d_status()` to `lib/validation.sh`**
   - Check for VCF analysis outputs
   - Return: "Complete", "Partial", "Not Started"
   - Check for expected file patterns (consult `modules/step1d/` to understand outputs)

3. **Update `analyze_pipeline_status()` in `lib/pipeline_common.sh`**
   - Add Step 1C status check
   - Add Step 1D status check
   - Update status determination logic to include all 4 steps
   - Update return values to include step1c_needed, step1d_needed cases

**Verification**:
- `check_step1c_status()` and `check_step1d_status()` work correctly
- `analyze_pipeline_status()` returns correct status for all 4 steps
- Main entry point correctly routes based on all step statuses

---

### **TASK 4: Update Documentation**

**Goal**: Ensure all documentation accurately reflects the improved pipeline.

**Steps**:

1. **Update `README.md`**
   - Add Steps 1C and 1D to overview
   - Update directory structure section to match target structure
   - Add CLI usage examples:
     ```bash
     # CLI usage
     bin/gatk_pipeline.sh -d dataset_name -s 1a
     bin/gatk_pipeline.sh -d dataset_name -s full
     
     # Interactive usage
     bin/gatk_pipeline.sh
     wrappers/interactive/step1a_interactive.sh
     ```
   - Document step dependencies (1A → 1B → 1C → 1D)
   - Update quick start section

2. **Create step-specific guides** (if they don't exist):
   - `docs/step1a_guide.md` - Step 1A detailed guide
   - `docs/step1b_guide.md` - Step 1B detailed guide
   - `docs/step1c_guide.md` - Step 1C detailed guide
   - `docs/step1d_guide.md` - Verify and update if exists

3. **Create architecture documentation**:
   - `docs/architecture.md` - Pipeline flow diagram and component relationships
   - `docs/development_guide.md` - How to add new steps or modify existing ones

**Verification**:
- README.md accurately describes all 4 steps
- Directory structure in README matches actual structure
- All examples in documentation work correctly

---

### **TASK 5: Configuration and Environment Setup**

**Goal**: Add robust environment validation and setup scripts.

**Steps**:

1. **Create `bin/setup_environment.sh`**
   - Check for conda environment (if used)
   - Verify required tools exist:
     - GATK (`gatk` command)
     - BWA (`bwa` command)
     - samtools (`samtools` command)
     - FastQC (`fastqc` command)
   - Validate paths from `config/environment.sh`:
     - Reference genome exists
     - Known sites VCF exists
     - Adapter file exists
   - Provide helpful error messages if anything is missing
   - Exit with appropriate codes

2. **Create `lib/config_validation.sh`**
   - Function `validate_pipeline_config()`:
     - Check all required paths exist
     - Validate SLURM account/partition are set
     - Warn about missing optional configs
     - Return success/failure status
   - Source this in main entry point before execution

3. **Update main entry point**
   - Call `validate_pipeline_config()` at startup
   - Provide clear error messages if validation fails
   - Suggest running `bin/setup_environment.sh` if issues found

**Verification**:
- `bin/setup_environment.sh` correctly identifies missing tools/configs
- `validate_pipeline_config()` catches common configuration errors
- Main entry point validates config before proceeding

---

### **TASK 6: Code Quality Improvements**

**Goal**: Ensure code quality and consistency across all scripts.

**Steps**:

1. **Run shellcheck on all scripts**
   - Fix all shellcheck warnings
   - Common issues to address:
     - SC1090, SC1091: Sourcing files (may need to exclude if paths are dynamic)
     - SC2086: Double quote variables
     - SC2155: Declare and assign separately

2. **Standardize script headers**
   - All scripts should have:
     - Shebang: `#!/bin/bash`
     - `set -euo pipefail` for error handling
     - Script description comment
     - Author/date if applicable

3. **Consistent error handling**
   - Use `error_exit()` function from libraries
   - Consistent error message format
   - Proper exit codes

**Verification**:
- `shellcheck` passes on all `.sh` files (or warnings are documented)
- All scripts follow consistent patterns
- Error handling is uniform across codebase

---

## ✅ VERIFICATION CHECKLIST

After completing all tasks, verify:

- [ ] `bin/gatk_pipeline.sh` exists and is executable
- [ ] `bin/gatk_pipeline.sh --help` shows usage information
- [ ] All 4 steps can be executed via CLI: `bin/gatk_pipeline.sh -d test -s {1a|1b|1c|1d}`
- [ ] Interactive mode works: `bin/gatk_pipeline.sh` (no args)
- [ ] All 4 interactive wrappers exist: `wrappers/interactive/step1{a,b,c,d}_interactive.sh`
- [ ] All interactive wrappers work standalone
- [ ] Status checking works for all 4 steps
- [ ] `archive/` directory contains redundant files
- [ ] No `.DS_Store` files in repository
- [ ] README.md documents all 4 steps
- [ ] Directory structure matches target structure
- [ ] All scripts pass shellcheck (or warnings documented)

---

## 📝 IMPORTANT NOTES FOR AI AGENTS

### **When Making Changes**

1. **Preserve existing functionality**
   - Don't break existing workflows
   - Test that old scripts still work (if kept for compatibility)
   - Maintain backward compatibility where possible

2. **Follow existing patterns**
   - Study `step1d_interactive.sh` for interactive script pattern
   - Study `modules/step1a/bin/run_step1a.sh` for module execution pattern
   - Follow logging patterns from `lib/logging.sh`

3. **Update path references**
   - When moving files, update all `SCRIPT_DIR` calculations
   - Use relative paths from script location, not hardcoded paths
   - Test path resolution works from different directories

4. **Error handling**
   - Always use `set -euo pipefail` in scripts
   - Provide clear, actionable error messages
   - Exit with appropriate codes (0 = success, non-zero = failure)

5. **Testing**
   - Test CLI flags work correctly
   - Test interactive mode still works
   - Test each step can be executed independently
   - Verify status checking returns correct values

### **Architecture Decisions**

1. **Wrapper naming**: Keep `wrappers/sbatch/` name (even though they execute, not submit directly) to avoid breaking changes. Document clearly in comments.

2. **Step integration**: Steps 1C and 1D should be optional in "full" pipeline run. User should be able to run 1A→1B→1C→1D sequentially, but each step should also work standalone.

3. **Interactive vs CLI**: Both interfaces should be available. CLI for automation, interactive for step-by-step guidance.

---

## 🚨 COMMON PITFALLS TO AVOID

1. **Hardcoded paths**: Always use `SCRIPT_DIR` and relative paths
2. **Breaking changes**: Don't remove functionality without deprecation period
3. **Inconsistent naming**: Follow existing conventions (step1a, step1b, etc.)
4. **Missing error handling**: Always validate inputs and handle errors gracefully
5. **Forgotten updates**: When moving files, update all references and documentation

---

**End of Implementation Guide**
