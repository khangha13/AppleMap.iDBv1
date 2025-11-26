# Apple GATK Modular Pipeline – Design, Incidents & AI Playbook

> **This file is now the single source of truth for `docs/`.  
> All other documents in `docs/` and `docs/step1d/` have been consolidated here.**

This document is written for **AI agents and future maintainers**. It captures:

- The **current architecture** of the modular GATK pipeline.
- The **authoritative path/variable model** (RDM, scratch, logs).
- The **logging + SLURM model**.
- The **behaviour and pitfalls of Steps 1A, 1B, 1C (VCF QC / Beagle)**.
- A curated **incident log with fixes** that explains why the design looks the way it does.

If you are an AI modifying this pipeline, assume this file is correct and that any
deleted `docs/*.md` files are superseded by what is written here.

---

## 1. High‑Level Architecture & Repo Layout

**User story:** a user wants to run a GATK pipeline (1A/1B/1C) from a *single* entry
point, while we keep the implementation modular and SLURM‑aware underneath.

- **Entry point:** `bin/gatk_pipeline.sh`
  - Reads config, resolves pipeline root, figures out which step(s) need to run.
  - Calls into `lib/` and `modules/step1[a–c]/` for actual work.
- **Libraries (`lib/`):**
  - `logging.sh` – shared professional logging primitives.
  - `slurm.sh` – helpers for job submission/status, `create_slurm_script`, etc.
  - `validation.sh` – input and completion validation helpers.
  - `pipeline_common.sh` – common pipeline utilities (RDM path resolution, etc.).
  - `env_bootstrap.sh` – pipeline root discovery and environment bootstrapping.
- **Step modules (`modules/`):**
  - `step1a/` – per‑sample variant calling (array over samples).
  - `step1b/` – `CombineGVCFs` across chromosomes into consolidated GVCFs.
  - `step1c/` – imputation / VCF QC tasks; see §4.3 and §4.4.
- **Wrapper scripts (`wrappers/`):**
  - Interactive / sbatch wrappers (e.g. `Interactive_1A.sh`, `Interactive_GATK.sh`,
    `step1b_submit.sh`, `step1c_submit.sh`), which build submission scripts and call
    into modules.

From a maintainer/AI perspective: **never replicate logic** that already exists in
`lib/` or `modules/`. Instead, add a thin wrapper and call shared helpers.

---

## 2. Authoritative Path & Variable Model

This section condenses everything from `RDM_PATH_CONSISTENCY_REPORT.md` and
`VARIABLE_CONSISTENCY_REPORT.md`. The TL;DR: **the pipeline’s path logic is now
consistent—do not re‑invent or partially duplicate it.**

### 2.1 RDM Base vs Dataset Paths (RDM invariants)

- **Environment (typically in `environment.template.sh`):**
  - `PIPELINE_RDM_BASE="/QRISdata/Q8367"`
  - `PIPELINE_RDM_DATASETS="${PIPELINE_RDM_BASE}/WGS_Reference_Panel"`
  - Shared resources (known sites, adapters, etc.) also hang off `PIPELINE_RDM_BASE`
    ❯ **never** from `PIPELINE_RDM_DATASETS`.
- **Config (`config/pipeline_config.sh`):**
  - `RDM_BASE_PATH="${PIPELINE_RDM_BASE:-...}"`
  - `RDM_DATASETS_PATH="${PIPELINE_RDM_DATASETS:-${RDM_BASE_PATH}/WGS_Reference_Panel}"`
- **Dataset inference (e.g. `lib/pipeline_common.sh`, `bin/gatk_pipeline.sh`):**
  - Dataset directory:
    - `${RDM_DATASETS_PATH}/${dataset_name}`  
      e.g. `/QRISdata/Q8367/WGS_Reference_Panel/NCBI`.
  - All step logic that needs a dataset path must *start from*
    `RDM_DATASETS_PATH`, **never** from `RDM_BASE_PATH` alone.

**Key invariant:** `RDM_BASE_PATH` describes the top‑level project root on RDM;
`RDM_DATASETS_PATH` describes the *panel root* (`WGS_Reference_Panel`); individual
dataset roots are under `RDM_DATASETS_PATH`.

### 2.2 Dataset path as a first‑class value (`DATASET_PATH`)

Template scripts for steps 1A/1B/1C originally accepted `$1` as
`RDM_BASE_PATH` even though they were actually handed a full dataset path. This
was confusing but functionally worked.

The templates have been standardised to:

- Parameter name: `DATASET_PATH`
- Value: full dataset path, e.g.  
  `/QRISdata/Q8367/WGS_Reference_Panel/NCBI`
- Usages inside templates:
  - `"${DATASET_PATH}/1.FASTQ"`
  - `"${DATASET_PATH}/5.Individual_VCF"`
  - `"${DATASET_PATH}/7.Consolidated_VCF"` (or module‑specific equivalents)

**AI rule:** when designing new templates or refactoring old ones, always
accept `DATASET_PATH` and treat it as the canonical root for per‑dataset
directories.

### 2.3 Scratch & log base paths

Derived from `PIPELINE_SCRATCH_BASE` and `PIPELINE_LOG_BASE` (in
`environment.template.sh` / user env):

- `SCRATCH_BASE_PATH` – root for temporary working directories on local scratch.
- `LOG_BASE_PATH` – root for pipeline and SLURM logs.

Defaults (can vary by site):

- `SCRATCH_BASE_PATH=/scratch/user/$USER`
- `LOG_BASE_PATH="${SCRATCH_BASE_PATH}/logs"`

**Invariants:**

- Every step should:
  - Use `SCRATCH_BASE_PATH` to create a working dir (module‑specific subdir).
  - Write logs under `${LOG_BASE_PATH}/${dataset_name}`.
- Avoid hard‑coding scratch under `/scratch/user/...` in modules; keep that in
  `PIPELINE_SCRATCH_BASE`/config.

---

## 3. SLURM & Logging Model

The logging/SLURM design is the synthesis of:

- `Logging_System_Guide.md`
- `SLURM_Aware_Logging_Summary.md`
- `Core_Pipeline_Logging_Guide.md`

### 3.1 Three‑tier logging

1. **Wrapper logs** – *interactive scripts* (`Interactive_1A.sh`, etc.)
   - **Location:** `./logs/wrapper/`
   - **Purpose:** user input, validation, dataset detection, step decisions,
     SLURM submission details.
   - Creates job tracking files like:
     - `./logs/wrapper/job_<JOB_ID>_<DATASET>.info` with:
       - dataset, sample count, SLURM script name,
       - recommended monitoring commands,
       - SLURM log directory.

2. **Pipeline (per‑sample / per‑step) logs** – *inside array jobs*
   - **Location:** `${LOG_BASE_PATH}/${dataset}/pipeline_<Sample>_*.log`
   - **Purpose:** structured, per‑sample step log:
     - start/end of each operation with durations,
     - SLURM context (job ID, array task),
     - fatal errors with context.
   - Uses logging functions like:
     - `pipeline_log_debug/info/warn/error/fatal`.
   - Legacy functions like `log_message` are shimmed to call the structured
     logging, so old code still writes into the unified format.

3. **SLURM job logs** – *stdout/stderr* of each job/array task
   - **Location:** `${LOG_BASE_PATH}/${dataset}/`
   - **Files:**
     - For Step 1A: `Apple_GATK_1A_*_*.output/.error`
     - For Step 1B: `1B_*_*_Apple_GATK_1B_*.output/.error` (or similar pattern).
   - **Purpose:** raw tool output and system errors (BWA, GATK, samtools, etc.).

**Design principle:** these three tiers are **complementary, not redundant**.
Never try to “mirror” SLURM logs in our logging system.

### 3.2 Log levels and configuration

For wrappers:

- `LOG_LEVEL` – controls verbosity (`DEBUG`, `INFO`, `WARN`, `ERROR`, `FATAL`).
- `LOG_TO_CONSOLE`, `LOG_TO_FILE` – allow disabling either channel.

For pipeline logs:

- `PIPELINE_LOG_LEVEL`
- `PIPELINE_LOG_TO_CONSOLE`
- `PIPELINE_LOG_TO_FILE`

AI guidance:

- Keep default levels at `INFO`.
- Use `DEBUG` for chatty, iterative development; do not ship new modules that
  only work at `DEBUG` without a good reason.

### 3.3 SLURM configuration consistency (Step‑specific variables)

From `VARIABLE_CONSISTENCY_REPORT.md`:

- **Steps 1A & 1B**:
  - Use step‑specific variables (e.g. `STEP1A_ACCOUNT`, `STEP1B_ACCOUNT`) that
    default to `PIPELINE_SLURM_ACCOUNT` but can be overridden.
  - Respect `SLURM_QOS` where provided.
- **Steps 1C & 1D (earlier state):**
  - Used `PIPELINE_SLURM_ACCOUNT` directly; did **not** respect QOS.
  - Config docs flag this as an inconsistency; recommended fix:
    - Add `STEP1C_ACCOUNT`, `STEP1D_ACCOUNT`, `STEP1C_PARTITION`,
      `STEP1D_PARTITION`, etc. to `environment.template.sh`.
    - Teach the `get_step1c_config` / `get_step1d_config` helpers to read these
      and optionally include `SLURM_QOS`.

**AI rule:** when touching Step 1C/1D config or adding new steps:

- Always match the Step 1A/1B pattern:
  - `STEPX_ACCOUNT="${STEPX_ACCOUNT:-$PIPELINE_SLURM_ACCOUNT}"`
  - `[ -n "${SLURM_QOS}" ] && echo "QOS=${SLURM_QOS}"`.

---

## 4. Step‑Specific Behaviour & Nuances

### 4.1 Step 1A – Per‑sample variant calling (`Interactive_1A.sh`)

Key points from `Interactive_1A_README.md`:

- **Role:** orchestrates per‑sample variant calling (originally
  `Apple_GATK_Pipeline_1a_CallVariantsPerSampleKHv2.sh`) via SLURM arrays.
- **Responsibilities:**
  - Prompt for dataset name and infer RDM dataset path:
    - `/QRISdata/Q8367/WGS_Reference_Panel/<dataset>`
  - Validate standard RDM directory structure:
    - `1.FASTQ/`, `2.FASTQC_pre_trimmed/`, `3.FASTQC_post_trimmed/`,
      `4.BAM/`, `5.Individual_VCF/`, `6.Consolidated_GVCF/`.
  - Generate or accept a sample list.
  - Perform a **completion check**:
    - Completed: has final per‑sample VCF + index.
    - Partial: has GVCF but missing genotyped VCF.
    - Incomplete: no outputs.
  - Derive a SLURM array range `0..N-1` and build a submission script.
  - Log all decisions, completion summaries, and job IDs via wrapper logging.

**AI‑relevant invariants:**

- Directory names under the dataset root are **part of the contract**. Code that
  validates or constructs paths must align with the RDM tree described in
  `Data_management/RDM_Database_Tree` and reflected here.
- Completion logic is path‑driven; avoid embedding “magic” file naming rules in
  multiple places—centralise in completion checks.

### 4.2 Step 1B – CombineGVCFs automation

From `Step1B_Automation_Guide.md` and the modular design:

- **Role:** combines per‑sample GVCFs (from `5.Individual_VCF`) into
  per‑chromosome consolidated GVCFs under `6.Consolidated_GVCF/`.
- **Automation in `Interactive_GATK.sh`:**
  - Detects Step 1A and Step 1B completion status for a dataset.
  - Offers to run Step 1B automatically when 1A is complete but 1B missing.
  - Can generate a “monitor Step 1A, then run 1B” helper script that:
    - Polls the Step 1A job via SLURM.
    - On success, runs `run_step1b`.
    - Propagates exit codes so that failures surface properly.
- **Completion checks:**
  - Look for consolidated chromosome VCFs + `.tbi` under `6.Consolidated_GVCF/`.

**AI‑relevant invariants:**

- Step 1B is **chromosome‑array based**; array configuration (e.g. `0–17`)
  should always come from a single source of truth (config + reference FASTA’s
  `.fai` or explicit chromosome list).
- Logging and SLURM patterns should mirror Step 1A for consistency.

### 4.3 Step 1C – VCF QC & Analysis (1c_VCF_QC pipeline)

The `docs/step1d/` files (`VCF_Analysis_README.md`, `*_SETUP_*`, `*_WORKFLOW*`,
`HPC_*`) all describe a **VCF QC pipeline** that is integrated as “Step 1C”
conceptually:

- **Inputs:** chromosome‑wise VCFs, e.g. `Chr00.vcf.gz` … `Chr17.vcf.gz`.
- **High‑level workflow:**

  1. Extract per‑site mean depth across samples into `SNP_site_meanDP.tsv`.
  2. Generate per‑chromosome depth PDFs (coverage vs position).
  3. Generate per‑site missingness‑vs‑depth PNGs for each chromosome.
  4. Combine all missingness plots into a single grid PNG.

- **Core scripts (in `1c_VCF_QC` project):**
  - `master_vcf_analysis_array.sh` – orchestrates a 4‑stage job chain:
    - TSV extraction → depth plots → missingness (array) → combine.
  - `master_vcf_analysis.sh` – sequential, single‑job variant.
  - `extract_site_mean_DP.sh` – standalone TSV extraction.
  - `vcf_analysis_config.sh` – configuration file.
  - R scripts:
    - `Average_depth_per_Chromosome.R`
    - `per_site_missingness_vs_depth.R`
    - `combine_per_site_missingness_plots.R`

#### 4.3.1 HPC‑ready configuration model

From `HPC_SETUP_GUIDE.md`, `HPC_MIGRATION_SUMMARY.md`,
`HPC_DEPLOYMENT_CHECKLIST.md`, and `VCF_Analysis_SETUP_SUMMARY.md`:

- All **hard‑coded local paths** and user names were purged.
- `vcf_analysis_config.sh` exports a set of **overridable env vars**:
  - `VCF_DIR` – directory containing input VCFs.
  - `VCF_PATTERN` – format string for filenames, e.g. `"Chr%02d.vcf.gz"`.
  - `CHR_START`, `CHR_END` – chromosome index range (usually `0..17`).
  - `WORK_DIR` – where outputs go (often `${VCF_DIR}`).
  - `LOG_DIR` – where pipeline logs go.
  - `R_SCRIPTS_DIR` – usually resolved relatively via `SCRIPT_DIR/../Rscript`.
  - `SLURM_ACCOUNT`, `SLURM_PARTITION` – cluster‑specific.
  - `PARALLEL_JOBS` and `*_RESOURCES` – CPU/memory/time triplets.
  - `CONDA_ENV`, `MODULE_TO_LOAD` – cluster R environment module names.
- `validate_config` and `show_config` encapsulate sanity checks; all scripts
  source `vcf_analysis_config.sh` rather than re‑deriving paths.

**AI‑relevant invariants:**

- Treat `vcf_analysis_config.sh` as the **sole source of truth** for Step 1C
  (VCF QC) inputs/outputs and resource settings.
- For portability, never add new hard‑coded paths or user names—add env vars
  with defaults and document them here instead.

#### 4.3.2 Execution modes

- **Array pipeline (recommended on HPC):** `master_vcf_analysis_array.sh`
  - Stages:
    1. Extract TSV (`SNP_site_meanDP.tsv`).
    2. Run `Average_depth_per_Chromosome.R`.
    3. Run `per_site_missingness_vs_depth.R` for each chromosome as one array.
    4. Run `combine_per_site_missingness_plots.R`.
  - Benefits: parallel computation, lower wall‑clock time.
- **Single‑job pipeline:** `master_vcf_analysis.sh`
  - Useful for testing/small datasets; slower for large sets.
- **Manual/stepwise:** run `extract_site_mean_DP.sh`, then R scripts individually.

From an AI standpoint, the important point is **how the data flows** (VCF → TSV
→ depth PDFs + missingness PNGs → combined PNG) and where to plug in new QC or
visualisation steps if desired.

### 4.4 Step 1C – Imputation scratch & logs (Beagle job)

`STEP1C_SUBMISSION_GUIDE.md` describes how the imputation/Beagle Step 1C job
manages working directories and logs:

- The job template (`modules/step1c/templates/step1c_job.sh`) chooses its temp
  working dir as:

  ```bash
  WORK_TMPDIR="${TMPDIR:-$(mktemp -d "${SCRATCH_BASE_PATH%/}/step1c_XXXXXX")}"
  ```

  - If `TMPDIR` is set, it is used **directly**.
  - Otherwise, a new temp directory is created under `SCRATCH_BASE_PATH`.

- Logs:
  - **SLURM logs:** `${LOG_BASE_PATH}/${dataset}/1C_*.output/.error`.
  - **Pipeline logs:** `${LOG_BASE_PATH}/${dataset}/pipeline_step1c_*.log`.

**AI‑relevant invariants:**

- Respect the precedence:
  - `TMPDIR` > `SCRATCH_BASE_PATH` when deciding where to put heavy temporary
    working files.
- For reproducible behaviour across clusters, fold new scratch behaviour into
  this pattern instead of ad‑hoc `cd /scratch/...` calls in modules.

---

## 5. Major Incidents & Resolutions (Postmortem Log)

This section merges the original `past_problems_and_resolutions.md` with the
variable/path consistency reports and key Step 1C/1D design decisions.

### 5.1 Pipeline root resolution failures

**Symptoms**

- Jobs running under `/var/spool/slurmd/...` failed with:
  - `Unable to locate pipeline root`
  - `No such file or directory` when sourcing `lib/logging.sh` or friends.
- Repo may live under `$HOME`, scratch, or RDM; compute nodes are only
  guaranteed `$HOME`.
- Some sbatch invocations did not include `PIPELINE_ROOT`, so child scripts
  guessed paths from their own working directory (the spool dir).

**Fixes**

1. **`lib/env_bootstrap.sh`:**
   - Allows admins/users to define `PIPELINE_ROOT` and other key vars in
     `$HOME/.gatk_pipeline_env` (overridable via `PIPELINE_ENV_FILE`).
2. **`resolve_pipeline_root` search order in `bin/gatk_pipeline.sh`:**
   1. Explicit `PIPELINE_ROOT` (from environment).
   2. `SLURM_SUBMIT_DIR` and its parent.
   3. Script‑relative paths.
   4. `$PIPELINE_HOME_CANDIDATE` or
      `$HOME/${PIPELINE_DIR_NAME:-GATK_Pipeline_KH_v1}`.
3. **All SLURM submission helpers** (including `submit_self`, module
   `submit_job`, and wrappers) now export:
   - `--export=ALL,PIPELINE_ROOT=...` plus any helper variables.
4. All persistent documentation about pipeline root resolution was consolidated
   here.

**Current rule of thumb:**

- **Never** call `sbatch` for a pipeline component without explicitly exporting
  `PIPELINE_ROOT`.
- If you move the repo, regenerate any pre‑built SLURM scripts that hard‑code
  `${PIPELINE_ROOT}`.

### 5.2 Master loop waiting forever after failed submission (2025‑11‑21)

**Symptom**

- When a step submission (e.g. Step 1B) failed immediately, the “master” job
  continued polling job status, seeing only “Not Started” indefinitely.

**Root cause**

- `wait_until_complete` only looked at status functions, which reported
  “Not Started” when the underlying `sbatch` submission itself had already
  failed.

**Fix**

- Each `run_step1X` helper now:
  - Returns the wrapper’s **exit code** (propagating `sbatch` failures).
- The master loop:
  - Aborts if submission fails.
  - Stops polling when a status function reports “Failed”.

**AI guidance:** when adding new orchestrated steps, ensure:

- Submission helpers propagate their own exit codes.
- The master orchestrator treats non‑zero returns as fatal and does not enter
  indefinite polling loops.

### 5.3 Step 1B RDM path regression (2025‑11‑21)

**Symptom**

- `run_step1b.sh` logged:
  - `RDM Base Path: /QRISdata/Q8367`
- But the user had selected:
  - `/QRISdata/Q8367/WGS_Reference_Panel/<dataset>`
- Step 1B then failed validation because it expected
  `5.Individual_VCF` under the full dataset root.

**Root cause**

- `wrappers/sbatch/step1b_submit.sh`:
  - Stored the full dataset path in `RDM_BASE_PATH`.
  - Then sourced `config/pipeline_config.sh`, where `RDM_BASE_PATH` was
    re‑defined to the cluster‑wide base `/QRISdata/Q8367`.
  - The sbatch command therefore forwarded the truncated path.

**Fix**

- The user‑provided full dataset path is preserved in `DATASET_RDM_PATH` before
  sourcing the config.
- The wrapper:
  - Logs `DATASET_RDM_PATH` explicitly.
  - Passes `DATASET_RDM_PATH` to `run_step1b.sh` (or other module entrypoints).
- Step 1A was unaffected because its wrapper `exec`’d the module without
  re‑sourcing the config in a way that clobbered the dataset path.

**AI guidance:** any future refactoring must:

- Keep `RDM_BASE_PATH` for *site‑wide* roots and
  `RDM_DATASETS_PATH` / `DATASET_PATH` / `DATASET_RDM_PATH` for dataset roots.
- Avoid reusing `RDM_BASE_PATH` to mean “this particular dataset’s root”.

### 5.4 Step 1B “No chromosomes detected” despite existing GVCFs (2025‑11‑21)

**Symptom**

- Step 1B aborted with:
  - `No chromosomes detected`
- But validated GVCFs clearly existed.
- The underlying problem: reference FASTA index missing on the compute node.

**Fix**

- Step 1B now ensures a **shared reference folder**, matching Step 1A, exists
  under:

  ```bash
  ${SCRATCH_BASE_PATH}/${dataset}_shared/Reference_genome
  ```

- On each run:
  - If this folder is missing, the orchestrator rsyncs `PIPELINE_REFERENCE_DIR`
    into place.
  - It ensures the staged FASTA exists.
  - If `.fai` is missing, it loads `samtools`, runs `samtools faidx`, and copies
    the index back to RDM for reuse.
  - Chromosome lists fall back to the staged `.fai` when
    `PIPELINE_CHROMOSOME_LIST` is empty.
- Failures show up cleanly via the master loop (see §5.2).

**AI guidance:**

- Any code that relies on chromosome enumeration should:
  - Prefer `.fai` and/or a single, central `PIPELINE_CHROMOSOME_LIST`.
  - Treat missing reference files/indexes as **configuration errors**, not as
    “no data”.

### 5.5 Step 1B SLURM script using `/var/spool` as repo root (2025‑11‑24)

**Symptom**

- Generated `step1b_array` scripts executed from `/var/spool/slurmd/...`.
- They tried to source `lib/functions.sh` relative to the current directory and
  failed with:
  - `No such file or directory`.

**Root cause**

- The generated SLURM script assumed the working directory at execution time
  would still be the repo root, which is not true once SLURM copies scripts to
  its spool.

**Fix**

- The script template now receives an absolute `PIPELINE_ROOT` placeholder.
- `create_slurm_script` substitutes the actual value during generation.
- The generated script sources libraries via:
  - `${PIPELINE_ROOT}/lib/...`

**AI guidance:** for any script that SLURM could copy into `/var/spool`:

- Always embed absolute paths (preferably via `${PIPELINE_ROOT}`).
- Do **not** rely on `$PWD` or relative paths to find `lib/` or `modules/`.

### 5.6 Variable and SLURM config inconsistencies (Steps 1C/1D)

From `VARIABLE_CONSISTENCY_REPORT.md`:

- **Issue:** Steps 1C and 1D used `PIPELINE_SLURM_ACCOUNT` directly and ignored
  QOS; Steps 1A/1B had step‑specific overrides and optional QOS support.
- **Impact:** Could not override accounts/partitions/QOS per step for 1C/1D.

**Recommended fix (may be partially or fully implemented depending on branch):**

- Add the following to `environment.template.sh` (or equivalent):
  - `STEP1C_ACCOUNT`, `STEP1C_PARTITION`, `STEP1C_NODES`, `STEP1C_NTASKS`.
  - `STEP1D_ACCOUNT`, `STEP1D_PARTITION`, `STEP1D_NODES`, `STEP1D_NTASKS`.
- Update `get_step1c_config` / `get_step1d_config` to:
  - Use these step‑specific variables, defaulting to `PIPELINE_SLURM_ACCOUNT`.
  - Include `SLURM_QOS` if set.

**AI guidance:**

- If you find any residual direct references to `PIPELINE_SLURM_ACCOUNT` inside
  step‑specific config helpers, consider them technical debt and align them with
  the Step 1A/1B pattern when safe to do so.

### 5.7 Step 1C/VCF QC HPC migration decisions

From `HPC_MIGRATION_SUMMARY.md` and companions:

- **What changed:**
  - All macOS‑specific and personal paths (e.g. `/Users/khangha/...`,
    `uqpha1`) were removed from scripts and kept only as examples in docs.
  - Config converted to parameterised env vars with sensible defaults.
  - R scripts referenced via relative paths using `SCRIPT_DIR` detection.
- **Why it matters:**
  - The VCF QC / Step 1C pipeline is now portable across HPC clusters.
  - Every new cluster needs only `vcf_analysis_config.sh` edits, not code
    changes.

**AI guidance:** when making future HPC‑targeted changes:

- Prefer config + validation functions over embedding cluster‑specific logic in
  bash or R.
- Only update this file to describe new required variables; do **not** reintroduce
  cluster‑specific paths into core scripts.

---

## 6. Configuration & Troubleshooting Checklist

This is the condensed, pipeline‑wide checklist that replaces the separate
`*_CHECKLIST`, `*_GUIDE`, and `*_SUMMARY` docs.

### 6.1 Core pipeline configuration (GATK 1A/1B/1C)

- **Per‑user/site environment:**
  - `~/.gatk_pipeline_env` (or `PIPELINE_ENV_FILE`) should export at least:
    - `PIPELINE_ROOT=/abs/path/to/GATK_Pipeline_KH_v1`
    - `PIPELINE_RDM_BASE=/QRISdata/Q8367`
    - `PIPELINE_RDM_DATASETS=/QRISdata/Q8367/WGS_Reference_Panel`
    - `PIPELINE_SCRATCH_BASE=/scratch/user/$USER`
    - `PIPELINE_LOG_BASE=/scratch/user/$USER/logs` (optional override).
- **Config sanity:**
  - `RDM_BASE_PATH` derived from `PIPELINE_RDM_BASE`.
  - `RDM_DATASETS_PATH` derived from `PIPELINE_RDM_DATASETS` with fallback.
  - All dataset paths built from `RDM_DATASETS_PATH` + dataset name.
- **SLURM:**
  - `PIPELINE_SLURM_ACCOUNT` set.
  - Step‑specific overrides optional but supported for 1A/1B (and recommended
    for 1C/1D once implemented).
  - `SLURM_QOS` respected where implemented (1A/1B; extend to others).

### 6.2 Submission wrappers

- Never call `sbatch` for pipeline jobs without:
  - `--export=ALL,PIPELINE_ROOT=...` (plus any step‑specific env).
- For new wrappers:
  - Use existing helpers in `lib/slurm.sh` rather than building raw `sbatch`
    strings.
  - Ensure logs end up under `${LOG_BASE_PATH}/${dataset}`.

### 6.3 Common troubleshooting patterns

- **Symptom:** `Unable to locate pipeline root` on compute node
  - Check `PIPELINE_ROOT` and that `$PIPELINE_ROOT/lib/logging.sh` exists.
  - Ensure `$PIPELINE_ENV_FILE` is accessible from compute nodes.
- **Symptom:** module array scripts launch with wrong repo path
  - Verify generated SLURM scripts contain a `PIPELINE_ROOT` line.
  - Remove stale scripts under `${PIPELINE_SLURM_SCRIPT_DIR}` and regenerate.
- **Symptom:** Step 1B says “No chromosomes detected”
  - Confirm the reference FASTA and `.fai` exist under
    `${SCRATCH_BASE_PATH}/${dataset}_shared/Reference_genome`.
  - Rebuild from `PIPELINE_REFERENCE_DIR` if needed.
- **Symptom:** Step 1C VCF QC jobs fail on HPC
  - Source `vcf_analysis_config.sh` and run `validate_config`.
  - Check that `VCF_DIR`, `VCF_PATTERN`, `LOG_DIR`, `R_SCRIPTS_DIR`,
    `CONDA_ENV`, `MODULE_TO_LOAD` are correct.

---

## 7. Future Work (Design‑Level, Not Tasks)

These items summarise the “future work” sections from the consolidated docs.
They are **not** promises, just design ideas worth considering:

- Extend `env_bootstrap.sh` to preload frequently overridden paths (scratch,
  logs) and to generate a skeleton `~/.gatk_pipeline_env`.
- Add CI tests that:
  - Run `bin/gatk_pipeline.sh` from `/tmp` with only `$HOME` available.
  - Simulate missing reference indexes and ensure user‑friendly failures.
- Bring Steps 1C/1D fully in line with the Step 1A/1B SLURM config pattern
  (step‑specific accounts/partitions, QOS handling).
- Add an interactive configuration wizard for the VCF QC (Step 1C) pipeline
  that writes a validated `vcf_analysis_config.sh`.
- Add automated resource‑usage reporting (e.g. ingest `seff` output) to
  recommend better default `*_RESOURCES` values.

If you are an AI agent planning changes, prefer to *extend* this document rather
than creating new standalone docs—this file is the canonical context for
pipeline behaviour and historical decisions.
