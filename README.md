# GATK Pipeline KH v1

Production-ready, modular GATK variant calling for apple genomes with Slurm-aware logging and a single entry point.

## Overview
- Implements GATK Best Practices: per-sample calling → joint genotyping → Beagle imputation → QC/PCA.
- Single launcher `bash bin/gatk_pipeline.sh` supports auto-mode, specific steps, interactive prompts, and self-submission to Slurm.
- Status-aware: inspects the dataset layout (1.FASTQ … 9.Imputation_QC) and only runs missing steps; Step 1A can resume or target incomplete samples and rebuild filtered sample lists.
- Slurm resources, reference paths, and directories are centralized in `config/pipeline_config.sh` (fed by `config/environment.sh`).
- Master runs write per-run log bundles under `${LOG_BASE_PATH}/${dataset}_<timestamp>`; Step 1B also tracks running/failed flags for monitoring.
- Step 1D QC/PCA templates load `miniforge/25.3.0-3`, `bcftools`, and `plink/2.00a3.6-gcc-11.3.0` automatically.

## Quick Start
1) Configure the environment
```bash
cp config/environment.template.sh config/environment.sh
nano config/environment.sh   # set PIPELINE_RDM_BASE, PIPELINE_SCRATCH_BASE, references, SLURM defaults
```
`pipeline_config.sh` will fall back to the template (with a warning) if `environment.sh` is missing—copying it is required for production.

2) Launch the pipeline (recommended self-submit)
```bash
# Auto-detect what to run, submit the master controller to Slurm, and poll between steps
bash bin/gatk_pipeline.sh -d <dataset> --submit

# Force interactive prompts (prints a status report and asks which steps to run; submits if needed)
bash bin/gatk_pipeline.sh -i

# Run a specific step without the interactive helper
bash bin/gatk_pipeline.sh -d <dataset> -s 1b      # 1a/1b/1c/1d/full/auto
```
- Override the dataset root with `--rdm-path /custom/path`.
- Preview actions without submitting jobs using `--dry-run`.
- Master poll interval is `PIPELINE_MASTER_POLL_SECS` (defaults to `MONITOR_INTERVAL` from `pipeline_config.sh`, currently 30s).

## Data Layout (RDM)
```
/QRISdata/Q8367/WGS_Reference_Panel/{dataset}/
├── 1.FASTQ/
├── 2.FASTQC_pre_trimmed/
├── 3.FASTQC_post_trimmed/
├── 4.BAM/
├── 5.Individual_VCF/
├── 6.genomicsdb_output/
├── 7.Consolidated_VCF/
├── 8.Imputated_VCF_BEAGLE/
└── 9.Imputation_QC/
```
`RDM_DATASETS_PATH` defaults to `${PIPELINE_RDM_BASE}/WGS_Reference_Panel` (set in `config/environment.sh`).

## Pipeline Steps
### Step 1A – Per-sample calling
- FastQC → Trimmomatic → FastQC → BWA → MarkDuplicates → BQSR → HaplotypeCaller → GenotypeGVCFs.
- Creates backups after key stages (default steps 3/4/5/6) in `${PIPELINE_WORK_DIR}/step1a/<dataset>` to support resume.
- Detects complete/partial samples; can restrict to incomplete samples.
- Submit: `bash wrappers/sbatch/step1a_submit.sh <dataset> <rdm_base> [--sample <id> | --sample-list <file>]`.
- Defaults (configurable): 10 CPUs, 32G, 200h, array limit 100.
- Outputs: per-sample GVCF (`*_raw.g.vcf.gz`) and genotyped VCF (`*_genotyped.vcf.gz` + index) in `5.Individual_VCF/`; FastQC HTML/ZIP in `2.*` and `3.*`; BAMs and indexes in `4.BAM/`.

### Step 1B – Joint genotyping
- Builds a sample map, runs `GenomicsDBImport` per chromosome (`PIPELINE_CHROMOSOME_LIST`, default Chr00–Chr17), then `GenotypeGVCFs`; outputs to `7.Consolidated_VCF/`.
- Uses `${LOG_BASE_PATH}/${dataset}` for `step1b_running.flag` / `step1b_failed.flag` when driven by the master script.
- Submit: `bash wrappers/sbatch/step1b_submit.sh <dataset> <rdm_base>`.
- Defaults: 6 CPUs, 36G, 200h, array limit 25.
- Outputs: per-chromosome consolidated VCFs + `.tbi` in `7.Consolidated_VCF/`.

### Step 1C – Beagle imputation
- Consumes consolidated VCFs from `7.Consolidated_VCF/`; writes imputed outputs to `8.Imputated_VCF_BEAGLE/` (override with `STEP1C_OUTPUT_DIR`).
- Submit: `bash wrappers/sbatch/step1c_submit.sh <dataset> <rdm_base>`.
- Defaults: 8 CPUs, 48G, 48h.
- Imputation toggle: `STEP1C_SELF_IMPUTE` (default `false`) controls Beagle's `impute=` flag (`false` = phasing-only; `true` = impute untyped markers).
- Outputs: Beagle-imputed VCFs + `.tbi` in `8.Imputated_VCF_BEAGLE/`.

### Step 1D – QC and PCA
- Runs QC plots/metrics and optional PCA on per-chromosome VCFs.
- Expects VCFs in the provided directory (pattern `Chr%02d.vcf.gz` by default). Set `VCF_PATTERN` for custom names or export `STEP1D_EXPECTED_CHROMS` to change chromosome count (Beagle mode defaults to 17).
- Flags: `--beagle` (imputed metrics), `--pca-only`, `--remove-relatives` (requires `--pca-only`), `--dry-run`.
- Submit: `bash wrappers/sbatch/step1d_submit.sh <dataset> <vcf_dir> [--beagle] [--pca-only] [--remove-relatives] [--dry-run]`.
- Outputs (in `VCF_DIR` by default): `variant_site_metrics.tsv`, depth/missingness/quality plots, optional PCA folder (`STEP1D_PCA_DIR`, default `pca_analysis`).

## Module Launch Examples (direct and interactive)
Use these when you want to drive a single module yourself rather than letting the master script orchestrate all steps.

### Step 1A (per-sample calling)
- Slurm wrapper:  
  `bash wrappers/sbatch/step1a_submit.sh <dataset> <rdm_base> [--sample <id> | --sample-list <file>]`
  - `--sample <id>`: run only one sample (expects `<id>_1.fastq.gz` / `<id>_2.fastq.gz`).
  - `--sample-list <file>`: one basename per line; runs only those samples.
  - No flags: builds a sample list from `1.FASTQ/`, detects complete samples, and can restrict to incomplete ones.
- Interactive helper: `bash wrappers/interactive/step1a_interactive.sh`
  - Prompts for dataset and sample selection; respects the same optional `--sample` / `--sample-list` flags.

### Step 1B (joint genotyping)
- Slurm wrapper:  
  `bash wrappers/sbatch/step1b_submit.sh <dataset> <rdm_base>`
  - Resources come from `config/pipeline_config.sh` (array jobs per chromosome). No extra flags.
- Interactive helper: `bash wrappers/interactive/step1b_interactive.sh`
  - Prompts for dataset and submits the Step 1B orchestrator.

### Step 1C (Beagle)
- Slurm wrapper:  
  `bash wrappers/sbatch/step1c_submit.sh <dataset> <rdm_base>`
  - Consumes `7.Consolidated_VCF/`; output goes to `8.Imputated_VCF_BEAGLE/` unless `STEP1C_OUTPUT_DIR` is set.
  - Imputation toggle: export `STEP1C_SELF_IMPUTE=true` to run Beagle with `impute=true`; default `false` (phasing-only/impute=false). Interactive wrapper prompts `Phasing-only mode? (y/n)`.
- Interactive helper: `bash wrappers/interactive/step1c_interactive.sh`
  - Prompts for dataset and launches imputation.

### Step 1D (QC/PCA)
- Slurm wrapper:  
  `bash wrappers/sbatch/step1d_submit.sh <dataset> <vcf_dir> [--beagle] [--pca-only] [--remove-relatives] [--dry-run]`
  - `--beagle`: expect Beagle-imputed inputs; uses AF/DR2/IMP metrics and adjusts plots.
  - `--pca-only`: skip QC plots/metrics, run only PCA (output in `${STEP1D_PCA_DIR}`).
  - `--remove-relatives`: drop KING-related samples before PCA (requires `--pca-only`).
  - `--dry-run`: print planned actions; do not create files.
- Interactive helper:  
  `bash wrappers/interactive/step1d_interactive.sh [--dir=PATH] [--vcf=Chr00,Chr01,...] [--beagle] [--pca-only] [--remove-relatives]`
  - `--dir=PATH`: preselect the VCF directory (otherwise prompted).
  - `--vcf=LIST`: comma-separated basenames to include; otherwise auto-detects `Chr??` or `*_snps.vcf.gz`.
  - Auto-inferrs `VCF_PATTERN` when possible; exports `VCF_INCLUDE_FILENAMES` for the master script.
  - Respects `STEP1D_EXPECTED_CHROMS` to change chromosome count; defaults to 17 in Beagle mode, 18 otherwise.

## Interactive Wrappers
- `wrappers/interactive/step1a_interactive.sh` (optional `--sample=<id>` / `--sample-list=<file>`).
- `wrappers/interactive/step1b_interactive.sh`
- `wrappers/interactive/step1c_interactive.sh`
- `wrappers/interactive/step1d_interactive.sh [--dir=PATH] [--vcf=Chr00,Chr01,...] [--beagle] [--pca-only] [--remove-relatives]`
  - Auto-detects VCF basenames (including Beagle-filtered `*_snps.vcf.gz`), infers `VCF_PATTERN` when possible, and exports `VCF_INCLUDE_FILENAMES` for the master script.

## Logging & Monitoring
- Master log bundle: `${LOG_BASE_PATH}/${dataset}_<timestamp>/GATK_master_script_<jobid>.{output,error}` (created when using `--submit` or running inside Slurm).
- Step 1B tracking flags: `${LOG_BASE_PATH}/${dataset}/step1b_running.flag` and `step1b_failed.flag` (or under the master log dir when orchestrated by `bin/gatk_pipeline.sh`).
- Generated Slurm scripts live in `${PIPELINE_SLURM_SCRIPT_DIR}`; temporary work in `${PIPELINE_WORK_DIR}`.
- Check jobs: `squeue -u $USER`; cancel: `scancel <jobid>`.

## Configuration Reference (edit in `config/environment.sh`)
- Paths: `PIPELINE_ROOT`, `PIPELINE_RDM_BASE`, `PIPELINE_RDM_DATASETS`, `PIPELINE_SCRATCH_BASE`, `PIPELINE_LOG_BASE`, `PIPELINE_WORK_DIR`, `PIPELINE_SLURM_SCRIPT_DIR`.
- References: `PIPELINE_REFERENCE_DIR`, `PIPELINE_REFERENCE_FASTA`, `PIPELINE_KNOWN_SITES_VCF`, `PIPELINE_ADAPTER_FASTA`.
- Slurm defaults: `PIPELINE_SLURM_ACCOUNT`, `PIPELINE_SLURM_PARTITION`, optional `PIPELINE_SLURM_QOS`.
- Step resources (override as needed): `STEP1A_CPUS/MEMORY/TIME/ARRAY_LIMIT`, `STEP1B_CPUS/MEMORY/TIME/ARRAY_LIMIT`, `STEP1C_CPUS/MEMORY/TIME`, `STEP1D_CPUS/MEMORY/TIME`.
- Chromosomes: `PIPELINE_CHROMOSOME_LIST` (default `Chr00 Chr01 … Chr17`).
- Step 1D tools/options: `PLINK2_BIN`, `BCFTOOLS_BIN`, `STEP1D_PCA_DIR`, `STEP1D_PCA_SHOW_LABELS`, `STEP1D_PCA_LABEL_SIZE`, `STEP1D_PCA_USE_GGREPEL`, `STEP1D_PCA_ONLY`, `STEP1D_REMOVE_RELATIVES`.
- Master polling: `PIPELINE_MASTER_POLL_SECS` (falls back to `MONITOR_INTERVAL`, default 30s).

## Testing
Run after configuring `config/environment.sh`:
```bash
bash test/test_modules.sh
bash test/test_comprehensive.sh
```
Warnings about missing external tools during isolated function tests are expected.

## Support
- Consult `docs/` (logging, modular design, Slurm integration, Step 1B automation).
- Check master and step logs under `${LOG_BASE_PATH}`.
- Verify Slurm job state (`squeue`) and rerun with `--dry-run` to validate configuration.

## Acknowledgments
- Authors: Phu Khang Ha
- Contributors: Paulo Henrique da Silva, Lisa Pasipanodya, Shashi Goonetilleke, Daniel Edge-Garza, Elizabeth Ross
- Institution: QAAFI, UQ
- Based on GATK Best Practices for variant calling

Version: 1.0  
Last Updated: October 2025  
License: Research Use
