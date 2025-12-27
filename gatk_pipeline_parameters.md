GATK Pipeline KH v1 – Configuration & Parameters
================================================

This document summarises the key configuration variables and parameters that control the behaviour of **GATK Pipeline KH v1**.
It is meant as a quick reference for new users and for maintainers updating defaults across HPC environments.

The canonical source of configuration is:

- `config/pipeline_config.sh`
- Environment file (cluster-/project-specific): `config/environment.sh` (or `config/environment.template.sh`)

Where possible, the pipeline follows this precedence:

1. **Environment variable** (exported before running the pipeline)
2. **Value set earlier in the shell session**
3. **Default set in `pipeline_config.sh`**

If you want to change a default for all runs on a given cluster, change it in `pipeline_config.sh` (or in `environment.sh` if the variable is meant to be environment-specific).

---

Storage & Path Configuration
----------------------------

These variables define where data and temporary files are stored.

| Variable | Default / Description |
| --- | --- |
| `PIPELINE_ROOT` | Automatically set to the top-level directory of the pipeline checkout. Used to resolve all internal paths. |
| `RDM_BASE_PATH` | Default: `/QRISdata/PLEASE_SET` or value from `PIPELINE_RDM_BASE` / `RDM_BASE_PATH`. Base RDM (QRISdata) path containing all datasets. |
| `RDM_DATASETS_PATH` | Default: `${RDM_BASE_PATH}/WGS_Reference_Panel` or `PIPELINE_RDM_DATASETS`. Root under which datasets live (e.g. one directory per dataset). |
| `SCRATCH_BASE_PATH` | Default: `/scratch/user/$(whoami)` or `PIPELINE_SCRATCH_BASE`. Base scratch space used for working directories and SLURM scripts. |
| `LOG_BASE_PATH` | Default: `${SCRATCH_BASE_PATH}/logs` or `PIPELINE_LOG_BASE`. Root for pipeline log files (not SLURM stdout/stderr). |
| `WRAPPER_LOG_PATH` | Default: `${PIPELINE_ROOT}/logs/wrapper`. Location for interactive wrapper logs. |
| `PIPELINE_SLURM_SCRIPT_DIR` | Default: `${SCRATCH_BASE_PATH}/gatk_pipeline/slurm_scripts`. Where auto-generated SLURM scripts are written. |
| `PIPELINE_WORK_DIR` | Default: `${SCRATCH_BASE_PATH}/gatk_pipeline/work`. Default working directory for pipeline steps. |

Reference & Resource Files
--------------------------

These variables control the reference genome, known-sites VCF, and adapter sequences.

| Variable | Default / Resolution | Notes |
| --- | --- | --- |
| `REFERENCE_GENOME` | Priority: `PIPELINE_REFERENCE_FASTA` → existing `REFERENCE_GENOME` → `/path/to/reference.fasta` | Path to the primary reference FASTA used by BWA and GATK. |
| `PIPELINE_REFERENCE_DIR` | Not set by default | If set and `REFERENCE_GENOME` is still the placeholder, the pipeline infers `reference.fasta` inside this directory. |
| `REFERENCE_GENOME_INDEX` | `${REFERENCE_GENOME}.fai` | FASTA index file used by samtools and GATK. |
| `REFERENCE_GENOME_DICT` | `${REFERENCE_GENOME%.*}.dict` | GATK sequence dictionary for the reference. |
| `KNOWN_SITES_VCF` | Priority: `PIPELINE_KNOWN_SITES_VCF` → existing `KNOWN_SITES_VCF` → `/path/to/known_sites.vcf.gz` | Known-sites VCF for BQSR. Must be indexed (`.tbi` or `.idx`). |
| `KNOWN_SITES_INDEX` | `${KNOWN_SITES_VCF}.tbi` | Tabix index for the known-sites VCF (created automatically if missing). |
| `ADAPTER_FILE` | Priority: `PIPELINE_ADAPTER_FASTA` → existing `ADAPTER_FILE` → `/path/to/adapters/TruSeq2-PE.fa` | Adapter sequences for Trimmomatic. Copied/symlinked into working directories. |

SLURM Defaults (Global)
------------------------

Base SLURM configuration used by the pipeline and exposed to step-specific settings.

| Variable | Default | Description |
| --- | --- | --- |
| `PIPELINE_SLURM_ACCOUNT` | `a_qaafi_cas` | Default SLURM project/account. |
| `PIPELINE_SLURM_PARTITION` | `general` | Default partition/queue. |
| `PIPELINE_SLURM_QOS` | `""` | Optional QoS string (empty by default). |
| `PIPELINE_SLURM_NODES` | `1` | Default number of nodes. |
| `PIPELINE_SLURM_NTASKS` | `1` | Default number of tasks per job. |
| `SLURM_ACCOUNT` | Same as `PIPELINE_SLURM_ACCOUNT` | Compatibility alias. |
| `SLURM_PARTITION` | Same as `PIPELINE_SLURM_PARTITION` | Compatibility alias. |
| `SLURM_QOS` | Same as `PIPELINE_SLURM_QOS` | Compatibility alias. |
| `SLURM_NODES` | Same as `PIPELINE_SLURM_NODES` | Compatibility alias. |
| `SLURM_NTASKS` | Same as `PIPELINE_SLURM_NTASKS` | Compatibility alias. |

Step‑Specific SLURM Settings
----------------------------

Each major step can override CPUs, memory, time limits, and array settings.

### Step 1A – Per-sample variant calling

| Variable | Default | Description |
| --- | --- | --- |
| `STEP1A_CPUS_PER_TASK` | `10` | CPUs per SLURM task for Step 1A array jobs. |
| `STEP1A_MEMORY` | `32G` | Memory request per task. |
| `STEP1A_TIME_LIMIT` | `200:00:00` | Walltime limit for each array task. |
| `STEP1A_ARRAY_MAX` | `100` | Maximum array size (0–N‑1). |
| `STEP1A_ACCOUNT` | `PIPELINE_SLURM_ACCOUNT` | SLURM account for Step 1A. |
| `STEP1A_PARTITION` | `PIPELINE_SLURM_PARTITION` | Partition used for Step 1A. |
| `STEP1A_NODES` | `PIPELINE_SLURM_NODES` | Nodes per job. |
| `STEP1A_NTASKS` | `PIPELINE_SLURM_NTASKS` | Tasks per job. |

### Step 1B – Combine GVCFs / joint genotyping

| Variable | Default | Description |
| --- | --- | --- |
| `STEP1B_CPUS_PER_TASK` | `6` | CPUs per SLURM task for Step 1B array jobs. |
| `STEP1B_MEMORY` | `36G` | Memory per task. |
| `STEP1B_TIME_LIMIT` | `200:00:00` | Time limit per array task. |
| `STEP1B_ARRAY_MAX` | `25` | Maximum array size for chromosomes. |
| `STEP1B_ACCOUNT` | `PIPELINE_SLURM_ACCOUNT` | SLURM account. |
| `STEP1B_PARTITION` | `PIPELINE_SLURM_PARTITION` | Partition. |
| `STEP1B_NODES` | `PIPELINE_SLURM_NODES` | Nodes per job. |
| `STEP1B_NTASKS` | `PIPELINE_SLURM_NTASKS` | Tasks per job. |

### Step 1C – Beagle imputation (placeholders)

| Variable | Default | Description |
| --- | --- | --- |
| `STEP1C_CPUS_PER_TASK` | `8` | Suggested CPUs for Beagle jobs. |
| `STEP1C_MEMORY` | `48G` | Suggested memory for Beagle. |
| `STEP1C_TIME_LIMIT` | `48:00:00` | Time limit for imputation jobs. |
| `STEP1C_ACCOUNT` | `PIPELINE_SLURM_ACCOUNT` | SLURM account. |
| `STEP1C_PARTITION` | `PIPELINE_SLURM_PARTITION` | Partition. |
| `STEP1C_NODES` | `PIPELINE_SLURM_NODES` | Nodes per job. |
| `STEP1C_NTASKS` | `PIPELINE_SLURM_NTASKS` | Tasks per job. |
| `STEP1C_SELF_IMPUTE` | `false` | Controls Beagle `impute=` flag. `false` = phasing-only (`impute=false`); `true` = self-impute untyped markers (`impute=true`). |

Actual imputation logic and Beagle-specific tuning currently live in `modules/step1c/templates/step1c_job.sh`.

### Step 1D – VCF QC and PCA

| Variable | Default | Description |
| --- | --- | --- |
| `STEP1D_CPUS_PER_TASK` | `16` | CPUs for Step 1D (QC + PCA) jobs. |
| `STEP1D_MEMORY` | `256G` | Memory request (intended for heavy PCA/QC). |
| `STEP1D_TIME_LIMIT` | `72:00:00` | Walltime limit for Step 1D jobs. |
| `STEP1D_ACCOUNT` | `PIPELINE_SLURM_ACCOUNT` | SLURM account. |
| `STEP1D_PARTITION` | `PIPELINE_SLURM_PARTITION` | Partition. |
| `STEP1D_NODES` | `PIPELINE_SLURM_NODES` | Nodes per job. |
| `STEP1D_NTASKS` | `PIPELINE_SLURM_NTASKS` | Tasks per job. |

Chromosome Configuration
------------------------

These variables control how chromosomes are enumerated across the pipeline.

| Variable | Default | Description |
| --- | --- | --- |
| `PIPELINE_CHROMOSOME_LIST` | `Chr00 Chr01 … Chr17` | Space-separated list of chromosome identifiers used across the pipeline (18 apple chromosomes). |
| `PIPELINE_CHROMOSOMES` | Array derived from `PIPELINE_CHROMOSOME_LIST` | Convenience array used by Step 1B and other modules. |

Software Paths & Tool Selection
-------------------------------

These variables define which executables or JARs are used by each tool. They can be overridden by environment modules or explicit paths.

| Variable | Default | Description |
| --- | --- | --- |
| `GATK_COMMAND` | `gatk` | Main GATK 4 command (can be a wrapper script or module-provided binary). |
| `TRIMMOMATIC_JAR` | If `EBROOTTRIMMOMATIC` set: `${EBROOTTRIMMOMATIC}/trimmomatic-0.39.jar`, else `trimmomatic-0.39.jar` | Location of the Trimmomatic JAR file. |
| `BWA_PATH` | `bwa` | Path to BWA executable (`bwa mem`). |
| `SAMTOOLS_PATH` | `samtools` | Path to samtools executable. |
| `FASTQC_PATH` | `fastqc` | Path to FastQC executable. |
| `PLINK2_BIN` | `plink2` | Name/path of PLINK2 binary used for PCA and QC. |
| `BCFTOOLS_BIN` | `bcftools` | Name/path of `bcftools` binary. |

Note: many HPC environments provide these tools as modules.
The templates under `modules/step1a/templates/`, `step1b/templates/`, and `step1d/templates/` load specific tested module versions (e.g. `fastqc/0.11.9-java-11`, `gatk/4.3.0.0-gcccore-11.3.0-java-11`, `plink/2.00a3.6-gcc-11.3.0`).

Trimming & Alignment Parameters (Step 1A)
-----------------------------------------

These parameters control Trimmomatic and BWA‑MEM behaviour in Step 1A. They are intended to match the original Apple pipeline.

### Trimmomatic

| Variable | Default | Description |
| --- | --- | --- |
| `TRIMMOMATIC_LEADING` | `20` | Minimum quality required at the beginning of a read. Bases below this are trimmed. |
| `TRIMMOMATIC_TRAILING` | `20` | Minimum quality at the end of a read. Bases below this are trimmed. |
| `TRIMMOMATIC_SLIDINGWINDOW` | `3:15` | Sliding window size and required average quality (`window:avgQ`). |
| `TRIMMOMATIC_AVGQUAL` | `20` | Minimum average read quality. |
| `TRIMMOMATIC_MINLEN` | `36` | Minimum length of reads to keep after trimming. |
| `TRIMMOMATIC_ILLUMINACLIP_SETTINGS` | `2:30:3:1:True` | Parameters appended to `ILLUMINACLIP:<adapter_file>:...` controlling adapter clipping behaviour. |

These are used in `run_trimmomatic()` (`modules/step1a/lib/functions.sh`) via a command of the form:

```bash
java -Xmx<MEM> -jar "${TRIMMOMATIC_JAR}" PE \
  -threads <threads> \
  *1.fastq.gz *2.fastq.gz \
  <paired/unpaired outputs> \
  ILLUMINACLIP:<adapter_basename>:${TRIMMOMATIC_ILLUMINACLIP_SETTINGS} \
  LEADING:${TRIMMOMATIC_LEADING} \
  TRAILING:${TRIMMOMATIC_TRAILING} \
  SLIDINGWINDOW:${TRIMMOMATIC_SLIDINGWINDOW} \
  AVGQUAL:${TRIMMOMATIC_AVGQUAL} \
  MINLEN:${TRIMMOMATIC_MINLEN}
```

### BWA‑MEM

| Variable | Default | Description |
| --- | --- | --- |
| `BWA_MEM_FLAGS` | `-M` | Extra flags passed to `bwa mem` (default `-M` to mark shorter split hits as secondary). |

These flags are used by `run_bwa_alignment()`:

```bash
"${BWA_PATH}" mem ${BWA_MEM_FLAGS} -t <threads> -R <RG> <ref.fa> \
  <sample>_forward_paired.fastq.gz \
  <sample>_reverse_paired.fastq.gz \
  | "${SAMTOOLS_PATH}" sort -o "${sample}_sorted.bam"
```

GATK Parameters
---------------

### Batch/Parallelisation Settings

| Variable | Default | Description |
| --- | --- | --- |
| `GATK_BATCH_SIZE` | `50` | Batch size used by `GenomicsDBImport` (controls number of samples processed at once). |
| `GATK_READER_THREADS` | `4` | Number of reader threads for GATK commands that support it. |

These are used by Step 1B functions (`modules/step1b/lib/functions.sh`) for `GenomicsDBImport` and related tools.

### Core GATK Tools (Usage Summary)

The following tools are invoked through `GATK_COMMAND` (see step‑specific modules for exact CLI options):

- **MarkDuplicates** – removes duplicate reads and writes metrics.
- **BaseRecalibrator** – builds BQSR recalibration tables using `KNOWN_SITES_VCF`.
- **ApplyBQSR** – applies recalibration to BAM files.
- **HaplotypeCaller** – emits per-sample GVCFs (`-ERC GVCF`) with thread count `--native-pair-hmm-threads`.
- **GenotypeGVCFs** – genotypes:
  - Per-sample GVCFs (Step 1A final VCF)  
  - GenomicsDB workspaces (Step 1B consolidated VCFs).
- **IndexFeatureFile** – ensures known-sites VCFs are indexed (`.tbi` / `.idx`).

Logging Configuration
---------------------

These variables control logging verbosity and destination across the pipeline.

| Variable | Default | Description |
| --- | --- | --- |
| `DEFAULT_LOG_LEVEL` | `INFO` | Baseline log level (DEBUG, INFO, WARN, ERROR). |
| `PIPELINE_LOG_LEVEL` | `${DEFAULT_LOG_LEVEL}` | Effective log level for pipeline logs. |
| `LOG_TO_CONSOLE` | `true` | Whether to print logs to stdout/stderr. |
| `LOG_TO_FILE` | `true` | Whether to log to a file on disk. |
| `PIPELINE_VERBOSE_LOGGING` | `false` | If `true`, `pipeline_config.sh` prints additional diagnostics when sourced. |
| `PIPELINE_GLOBAL_DRY_RUN` | `false` | Global dry‑run toggle; intended to be honoured by high‑level wrappers. |

Backup & Shared Reference Configuration
---------------------------------------

These settings control automatic backup of intermediate results and shared reference management for array jobs.

| Variable | Default | Description |
| --- | --- | --- |
| `PIPELINE_ENABLE_BACKUP` | `true` | Enables backup behaviour in Step 1A helper functions. |
| `PIPELINE_SHARED_REF_TIMEOUT` | `120` | Seconds to wait for shared reference files to be initialised in array jobs. |
| `PIPELINE_SHARED_REF_POLL_INTERVAL` | `3` | Polling interval (seconds) when waiting for shared reference files. |
| `PIPELINE_BACKUP_STEPS` | `3,4,5,6` | Comma-separated list of step numbers eligible for backups. |

Validation Configuration
------------------------

Used by validation helpers to check data layout and required files before running the pipeline.

| Variable | Default / Contents | Description |
| --- | --- | --- |
| `REQUIRED_DIRS` | `("1.FASTQ" "2.FASTQC_pre_trimmed" "3.FASTQC_post_trimmed" "4.BAM" "5.Individual_VCF" "6.genomicsdb_output" "7.Consolidated_VCF" "8.Imputated_VCF_BEAGLE" "9.Imputation_QC")` | Directories expected under each dataset path. |
| `FASTQ_EXTENSIONS` | `("_1.fastq.gz" "_2.fastq.gz")` | Allowed FASTQ filename suffixes. |
| `GVCF_EXTENSIONS` | `("_raw.g.vcf.gz" "_raw.g.vcf.gz.tbi")` | GVCF + index suffixes. |
| `VCF_EXTENSIONS` | `("_genotyped.vcf.gz" "_genotyped.vcf.gz.tbi")` | Per-sample VCF + index suffixes. |
| `CONSOLIDATED_VCF_EXTENSIONS` | `(".vcf.gz" ".vcf.gz.tbi")` | Consolidated VCF + index suffixes. |

Monitoring & Notification Configuration
---------------------------------------

These variables control how often running jobs are polled and whether email notifications are sent.

| Variable | Default | Description |
| --- | --- | --- |
| `MONITOR_INTERVAL` | `30` | Seconds between monitoring polls for active jobs. |
| `MAX_MONITOR_ATTEMPTS` | `100` | Maximum number of monitoring iterations before giving up. |
| `EMAIL_NOTIFICATIONS` | `false` | Whether to send email notifications on important events. |
| `EMAIL_ADDRESS` | `""` | Recipient email address when notifications are enabled. |

Advanced Behaviour Flags
------------------------

These toggles govern advanced behaviour like resume support, backups, cleanup, and debugging.

| Variable | Default | Description |
| --- | --- | --- |
| `ENABLE_RESUME` | `true` | Enables resume support using a `resume_step.txt` file. |
| `RESUME_STEP_FILE` | `resume_step.txt` | File that tracks the last completed step. |
| `ENABLE_BACKUPS` | `true` | Enables periodic backups of working directories. |
| `BACKUP_INTERVAL` | `5` | Backup interval in “logical steps” or pipeline units. |
| `CLEANUP_TEMP_FILES` | `true` | Whether temporary files are deleted after successful completion. |
| `CLEANUP_INTERMEDIATE_FILES` | `false` | Whether intermediate analysis files are removed. |
| `DEBUG_MODE` | `false` | Enables additional debug logging and verbose output. |
| `TEST_MODE` | `false` | Marks a run as test-mode (behaviour is script-dependent). |
| `VERBOSE_OUTPUT` | `false` | Enables extra user-facing output. |

Step 1D PCA & Plotting Options
------------------------------

Step 1D uses a combination of `bcftools`, `plink2`, and R scripts to generate QC plots and perform PCA.
Some of these options are currently configured directly in `modules/step1d/templates/master_vcf_analysis.sh` and `plink2_PCA.sh`.

The following variables already exist in `pipeline_config.sh`:

| Variable | Default | Description |
| --- | --- | --- |
| `STEP1D_PCA_ONLY` | `false` | If `true`, Step 1D runs PCA-only (QC stages may be skipped). |
| `STEP1D_REMOVE_RELATIVES` | `false` | If `true`, PCA helper can be asked to remove close relatives (configurable via CLI). |
| `STEP1D_PCA_DIR` | `pca_analysis` | Directory name inside the working directory where PCA outputs are written. |
| `STEP1D_PCA_SHOW_LABELS` | `true` | Default for whether sample labels are shown on PCA plots. |
| `STEP1D_PCA_LABEL_SIZE` | `1.5` | Default label size in PCA plots. |
| `STEP1D_PCA_USE_GGREPEL` | `true` | If `true`, PCA plots use `ggrepel` to reduce label overlap. |
| `STEP1D_PCA_MERGED_PATTERN` | `*merged*.vcf.gz,*merge*.vcf.gz` | Comma-separated glob(s) to locate a pre-merged multi-chromosome VCF for PCA (candidates containing `Chr` are ignored by default). |
| `STEP1D_PCA_FORCE_CONCAT` | `false` | If `true`, always concatenate per-chromosome VCFs and skip merged VCF detection. |
| `STEP1D_PCA_MERGED_EXCLUDE_CHR` | `true` | If `true`, ignore merged candidates containing `Chr` (case-insensitive) in the filename. |
| `STEP1D_DUPLICATE_MODE` | `flag` | KING-based duplicate detection mode: `off`, `flag` (report/highlight), or `remove` (exclude from PCA input). |
| `STEP1D_DUPLICATE_KING_THRESHOLD` | `0.45` | Kinship threshold for duplicate detection (pairs ≥ threshold are flagged). |
| `STEP1D_AF_PLOTS_DIR` | `af_distribution_plots` | Output directory name for allele frequency distribution plots. |
| `STEP1D_AF_HIST_BINS` | `50` | Histogram bin count for allele frequency plots. |

Script-level PCA and QC thresholds (not yet centralised) include:

- PLINK2 QC thresholds: `--geno 0.05`, `--mind 0.10`, `--maf 0.01`
- LD pruning: `--indep-pairwise 200 50 0.2`
- KING cutoff: `--king-cutoff 0.125`
- Duplicate detection threshold: `STEP1D_DUPLICATE_KING_THRESHOLD` (default `0.45`)
- Number of principal components: `--pca 20`
- Plot image format: `PLOT_IMAGE_FORMAT=png`
- Call-rate heatmap bins: `CALL_RATE_HEATMAP_BINS=100`

These are good candidates to be moved into `pipeline_config.sh` in a future refactor so that PCA/QC behaviour is fully configurable.

Helper Accessors
----------------

For use in scripts and templates, `pipeline_config.sh` provides helper functions that print out consolidated settings:

| Function | Description |
| --- | --- |
| `get_step1a_config` | Prints key Step 1A SLURM settings (CPUS, MEMORY, TIME, ACCOUNT, PARTITION, NODES, NTASKS, ARRAY_MAX, QOS). |
| `get_step1b_config` | Prints Step 1B SLURM settings in the same format. |
| `get_step1c_config` | Prints Step 1C SLURM settings. |
| `get_step1d_config` | Prints Step 1D SLURM settings. |
| `get_reference_fasta` | Echoes the effective reference FASTA path. |
| `get_known_sites_vcf` | Echoes the effective known-sites VCF path. |
| `get_adapter_fasta` | Echoes the effective adapter FASTA path. |
| `get_rdm_datasets_path` | Echoes the base datasets path on RDM. |
| `get_scratch_base_path` | Echoes the base scratch path. |
| `get_log_base_path` | Echoes the base log path. |
| `get_pipeline_env_source` | Shows which environment file (`environment.sh` or template) was sourced. |
| `pipeline_config_summary` | Prints a short multi-line summary (env source, RDM base, RDM datasets path, scratch base, reference FASTA). |

How to Override Parameters
--------------------------

There are three common ways to change pipeline behaviour:

1. **Environment file (`config/environment.sh`)**  
   - Best for cluster- or project-wide settings such as paths to reference genomes, RDM base paths, and default scratch directories.  
   - Example:
     ```bash
     export PIPELINE_REFERENCE_FASTA=/QRISdata/Q8367/Reference_Genome/GDDH13_1-1_formatted.fasta
     export PIPELINE_RDM_BASE=/QRISdata/Q8367
     export PIPELINE_SCRATCH_BASE=/scratch/user/$USER
     ```

2. **Environment variables at submission time**  
   - Best for one-off overrides (e.g. more RAM for a single dataset, changing GATK binary).  
   - Example:
     ```bash
     export STEP1A_MEMORY=48G
     export GATK_COMMAND=gatk-4.3.0.0
     sbatch wrappers/sbatch/step1a_submit.sh ...
     ```

3. **Editing `config/pipeline_config.sh`**  
   - Best for changing the **default behaviour** of the pipeline (e.g. adopting new QC thresholds or tool versions across the board).  
   - Changes here will affect all subsequent runs unless overridden by environment variables.

When in doubt, prefer setting values in `environment.sh` for site-specific paths and in `pipeline_config.sh` for logical defaults, and reserve per-job exports for exceptional cases.
