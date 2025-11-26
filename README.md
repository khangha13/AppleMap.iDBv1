# GATK Pipeline KH v1

A production-ready, modular GATK variant calling pipeline designed for apple genome analysis.

## 🧬 Overview

This pipeline implements the GATK Best Practices workflow for variant calling in apple genomes, featuring:

- **Step 1A**: Per-sample variant calling (HaplotypeCaller)
- **Step 1B**: Combine GVCFs (GenomicsDBImport + GenotypeGVCFs)
- **Step 1C**: Beagle imputation for missing genotypes
- **Step 1D**: VCF quality control and analysis
- **Professional logging system** with SLURM integration
- **Modular architecture** for easy maintenance and extension
- **Single entry point** with user-friendly interface (CLI and interactive modes)

## 🚀 Quick Start

**⚠️ Important: Before running the pipeline, you must set up your environment configuration (see [Configuration](#-configuration) section below).**

```bash
# 1. Copy and configure environment file
cp config/environment.template.sh config/environment.sh
nano config/environment.sh  # Edit with your specific paths and settings

# 2. Run the pipeline
bash bin/gatk_pipeline.sh -d <dataset_name>          # run locally
# Or submit to SLURM:
sbatch bin/gatk_pipeline.sh -d <dataset_name>
```

> ℹ️  All automation must invoke `bin/gatk_pipeline.sh`; the legacy entry point without the `.sh` suffix has been removed.

Each run produces a dedicated log folder named `<dataset>_<YYYYMMDD>` under your configured `LOG_BASE_PATH`. When you let the pipeline self-submit (`bin/gatk_pipeline.sh --submit` or via the interactive launcher), the master script redirects both stdout/stderr and the structured logger into that folder so every Slurm job has an isolated log bundle.

### **Interactive Mode Highlights**
- Running `bin/gatk_pipeline.sh` with no flags enters an interactive helper that now **prompts you about any pending steps (1B/1C/1D)** before it submits work to Slurm.
- You’ll see `Would you like to run Step 1B now? [y/n]` (etc.) directly in the login shell; answering `y` records the choice and launches the master job, while `n` exits without submission.
- If all expected outputs are already present, the interactive helper immediately reports “✅ No actions required” and exits—no Slurm submission is created.
- Once submitted, the master job honours the recorded action and continues through the remaining steps automatically, with status updates written to the per-run log directory noted above.

The pipeline will guide you through:
1. Dataset configuration
2. Pipeline status analysis
3. Automatic step detection and execution
4. Progress monitoring
5. (When run via sbatch) Automatic polling between steps until completion

## 🔧 Configuration

### **Initial Setup: Environment Configuration**

**⚠️ Required**: Before running the pipeline or tests, you must configure your environment:

```bash
# 1. Copy the template to create your environment file
cp config/environment.template.sh config/environment.sh

# 2. Edit environment.sh with your specific settings:
#    - PIPELINE_ROOT: Absolute path to the cloned repo (can live in $HOME, RDM, or scratch)
#    - PIPELINE_RDM_BASE: Base path for RDM storage
#    - PIPELINE_SCRATCH_BASE: Base path for scratch/logs
#    - PIPELINE_REFERENCE_DIR: Directory containing reference files
#    - PIPELINE_REFERENCE_FASTA: Path to reference genome FASTA
#    - PIPELINE_KNOWN_SITES_VCF: Path to known sites VCF for BQSR
#    - PIPELINE_ADAPTER_FASTA: Path to adapter sequences
#    - SLURM account, partition, and other SLURM defaults

nano config/environment.sh  # or use your preferred editor

# 3. Verify configuration by running tests (recommended)
bash test/test_modules.sh
bash test/test_comprehensive.sh
```

### **Environment Variables**
- `NUM_THREADS`: Number of CPU threads (default: 10)
- `MEMORY`: Memory allocation (default: 32G for Step 1A, 36G for Step 1B)
- `REFERENCE_GENOME`: Path to reference genome
- `KNOWN_SITES`: Path to known variant sites for BQSR
- `RDM_DATASETS_PATH`: Base path for datasets (set in `config/environment.sh`)
- `SCRATCH_BASE_PATH`: Scratch/log base (set in `config/environment.sh`)
- `PIPELINE_MASTER_POLL_SECS`: Master controller polling interval (default 600)

### **Run-Specific Log Folders**
- Set `PIPELINE_ROOT`, `PIPELINE_RDM_BASE`, and `PIPELINE_SCRATCH_BASE` to absolute paths that make sense for your site (e.g. `$HOME`, shared RDM, or scratch).
- `LOG_BASE_PATH` (defaults to `${SCRATCH_BASE_PATH}/logs`) is used to create per-run folders named `<dataset>_<YYYYMMDD>`; the master script streams stdout/stderr and the structured `logging.sh` output into that folder.
- You can relocate logs to RDM or scratch simply by pointing `PIPELINE_LOG_BASE` in `config/environment.sh` at the desired filesystem—no code changes required.
- When using interactive mode, the helper prints the exact log folder before submitting so you can `tail` the job output immediately.

### **SLURM Configuration**
- **Step 1A**: Array jobs for parallel sample processing
- **Step 1B**: Array jobs for parallel chromosome processing
- **Resource allocation**: Automatically configured based on dataset size

### **Step 1D PCA Controls**
- Enable PCA-only mode in Step 1D with `--pca-only` (interactive) or `export STEP1D_PCA_ONLY=true`.
- Remove relatives before PCA with `--remove-relatives` (requires PCA-only mode) or `export STEP1D_REMOVE_RELATIVES=true`.
- When calling the core script directly: `master_vcf_analysis.sh [--beagle] [--pca-only] [--remove-relatives]`.

## 🚀 Advanced Usage

### **Custom Function Development**
1. Add new functions to `modules/step1a/functions.sh` or `modules/step1b/functions.sh`
2. Update the corresponding `main.sh` to call new functions
3. Test with `bash test/test_comprehensive.sh`

### **Pipeline Extension**
1. Create new module directories in `modules/`
2. Follow the existing pattern: `main.sh`, `functions.sh`, `config.sh`
3. Update `bin/gatk_pipeline.sh` to include new steps


## 📁 Directory Structure

```
GATK_Pipeline_KH_v1/
├── archive/                         # Archived legacy scripts/configs
├── bin/
│   └── gatk_pipeline.sh           # Single entry point
├── modules/
│   ├── step1a/
│   │   ├── bin/run_step1a.sh
│   │   ├── lib/functions.sh
│   │   └── templates/step1a_array.sh
│   ├── step1b/
│   │   ├── bin/run_step1b.sh
│   │   ├── lib/functions.sh
│   │   └── templates/step1b_array.sh
│   ├── step1c/
│   │   ├── bin/run_step1c.sh
│   │   ├── lib/functions.sh
│   │   └── templates/step1c_job.sh
│   └── step1d/
│       ├── bin/run_step1d.sh
│       ├── lib/
│       └── templates/master_vcf_analysis.sh
├── wrappers/
│   ├── interactive/
│   │   ├── step1a_interactive.sh
│   │   ├── step1b_interactive.sh
│   │   ├── step1c_interactive.sh
│   │   └── step1d_interactive.sh
│   └── sbatch/
│       ├── step1a_submit.sh
│       ├── step1b_submit.sh
│       ├── step1c_submit.sh
│       └── step1d_submit.sh
├── config/
│   ├── environment.template.sh
│   └── pipeline_config.sh
├── lib/
│   ├── interactive_common.sh
│   ├── logging.sh
│   ├── pipeline_common.sh
│   ├── slurm.sh
│   └── validation.sh
├── docs/
├── test/
└── utils/
```

## 🔧 Features

### **Professional Logging System**
- Multi-level logging (DEBUG, INFO, WARN, ERROR, FATAL)
- SLURM-aware logging that complements SLURM's own output
- Separate log files for wrapper and pipeline execution
- Timestamped logs with professional formatting

### **Modular Architecture**
- **Step 1A Module**: Per-sample processing with 8 pipeline functions
- **Step 1B Module**: GVCF combination with 8 GATK/utility functions
- **Step 1C Module**: Beagle imputation and post-imputation preparation utilities
- **Step 1D Module**: QC, plotting, and optional PLINK2 PCA (with relatedness filtering)
- **Shared Libraries**: Reusable logging, SLURM, and validation functions
- **Easy Extension**: Add new steps or modify existing ones easily

### **User-Friendly Interface**
- Single command execution: `bash bin/gatk_pipeline.sh`
- Interactive dataset configuration
- Automatic pipeline status detection
- Intelligent completion checking
- Progress monitoring and job management

### **SLURM Integration**
- Dynamic SLURM script generation
- Array job support for parallel processing
- Automatic resource allocation
- Job monitoring and status tracking

## 📋 Pipeline Steps

### **Step 1A: Per-Sample Variant Calling**
1. **FastQC** on raw reads
2. **Trimmomatic** adapter removal and quality filtering
3. **FastQC** on trimmed reads
4. **BWA alignment** to reference genome
5. **GATK MarkDuplicates** for duplicate removal
6. **GATK BaseRecalibrator** for quality score recalibration
7. **GATK HaplotypeCaller** for variant calling (GVCF output)
8. **GATK GenotypeGVCFs** for VCF conversion

Note: After trimming (Step 2), the pipeline backs up only the paired trimmed FASTQs to scratch (unpaired reads are not backed up). This enables fast recovery without re-trimming if later steps fail.

### **Step 1B: Combine GVCFs**
1. **Archive creation** of individual GVCF files
2. **Sample map generation** for GenomicsDBImport
3. **GenomicsDBImport** per chromosome
4. **GenotypeGVCFs** for consolidated VCF generation
5. **Result copying** to permanent storage

## 🛠️ Usage Examples

### **Run Complete Pipeline**
```bash
# Interactive (prompts for dataset and routes automatically)
bash bin/gatk_pipeline.sh -i

# CLI auto mode for a specific dataset (runs next needed step)
bash bin/gatk_pipeline.sh -d NCBI_9sample_Noc2025

# Submit the master controller to SLURM (auto-run all steps w/ polling)
sbatch bin/gatk_pipeline.sh -d NCBI_9sample_Noc2025
```

### **Run Individual Steps**

#### **Using Direct SLURM Submission**
```bash
# Step 1A - Submit all samples as array job
bash wrappers/sbatch/step1a_submit.sh dataset_name /path/to/rdm/base

# Step 1A - Run a single specific sample
bash wrappers/sbatch/step1a_submit.sh dataset_name /path/to/rdm/base --sample YourSampleID

# Step 1A - Run a custom list (one basename per line)
bash wrappers/sbatch/step1a_submit.sh dataset_name /path/to/rdm/base --sample-list /abs/path/to/list.txt

# Step 1B - Submit chromosome array job
bash wrappers/sbatch/step1b_submit.sh dataset_name /path/to/rdm/base

# Step 1C - Submit imputation job
bash wrappers/sbatch/step1c_submit.sh dataset_name /path/to/rdm/base

# Step 1D - Submit QC and analysis job
bash wrappers/sbatch/step1d_submit.sh dataset_name /path/to/rdm/base
```

### **Interactive Launchers**

Interactive wrappers provide a user-friendly interface for running individual steps:

```bash
# Step 1A - Interactive mode (prompts for dataset and sample selection)
wrappers/interactive/step1a_interactive.sh
# Optional flags to preselect:
#   --sample=YourSampleID
#   --sample-list=/abs/path/to/list.txt

# Step 1B - Interactive mode
wrappers/interactive/step1b_interactive.sh

# Step 1C - Interactive mode
wrappers/interactive/step1c_interactive.sh

# Step 1D - Interactive mode
# Non-imputed consolidated VCFs (PCA-only)
wrappers/interactive/step1d_interactive.sh --dir="/.../7.Consolidated_VCF" --pca-only --remove-relatives

# Beagle-imputed VCFs (PCA-only)
wrappers/interactive/step1d_interactive.sh --dir="/.../8.Imputated_VCF_BEAGLE" --beagle --pca-only --remove-relatives
```

The Step 1D wrapper now auto-detects Beagle-filtered files such as `panel.snps.clean__Chr01_snps.vcf.gz`, infers the matching `VCF_PATTERN`, and proceeds without manual overrides. If your filenames cannot be inferred automatically, set a pattern explicitly (e.g. `export VCF_PATTERN='panel.snps.clean__Chr%02d_snps.vcf.gz'`) before launching.

When running with `--beagle`, the launcher automatically expects 17 chromosomes (Chr01–Chr17). If your imputed dataset uses a different chromosome count, override it via `export STEP1D_EXPECTED_CHROMS=<count>`.

During Step 1D execution the scripts automatically load `miniforge/25.3.0-3`, `bcftools`, and `plink/2.00a3.6-gcc-11.3.0`. Ensure those modules exist on your cluster (or pre-load equivalent modules/expose `PLINK2_BIN` before launching).

**Step 1D Flags:**
- `--pca-only`: bypass QC outputs and run only the PLINK2 PCA stage.
- `--remove-relatives`: remove close relatives before PCA (KING cutoff 0.125). Requires `--pca-only`.
- `--beagle`: treat inputs as Beagle-imputed (different metrics layout).

## 📊 Data Structure

The pipeline expects the following RDM directory structure:

```
/QRISdata/Q8367/WGS_Reference_panel/{dataset_name}/
├── 1.FASTQ/                     # Raw sequencing data
├── 2.FASTQC_pre_trimmed/        # Quality control (before trimming)
├── 3.FASTQC_post_trimmed/       # Quality control (after trimming)
├── 4.BAM/                       # Aligned and processed BAM files
├── 5.Individual_VCF/            # Per-sample GVCF files
├── 6.genomicsdb_output/         # GenomicsDB workspace (temporary, Step 1B)
├── 7.Consolidated_VCF/          # Combined VCF files (Step 1B output)
├── 8.Imputated_VCF_BEAGLE/      # Beagle-imputed VCF files (Step 1C output)
└── 9.Imputation_QC/             # Quality control metrics and plots (Step 1D output)
```

## 🔍 Monitoring and Logs

### **SLURM Job Monitoring**
```bash
# Check job status
squeue -u $USER

# Cancel job
scancel <job_id>

# View SLURM logs
ls /scratch/user/uqpha1/logs/{dataset_name}/
# Master controller default logs (if using SBATCH header in gatk_pipeline.sh):
#   /scratch/user/uqpha1/logs/GATK_master_script_<jobid>.output|error
```

### **Pipeline Logs**
- **Wrapper logs**: `./logs/wrapper/`
- **Pipeline logs**: `/scratch/user/uqpha1/logs/{dataset_name}/`
- **SLURM logs**: `/scratch/user/uqpha1/logs/{dataset_name}/`

## 📚 Documentation

- **Logging System Guide**: `docs/Logging_System_Guide.md`
- **Step 1B Automation**: `docs/Step1B_Automation_Guide.md`
- **Modular Implementation**: `docs/Modular_Implementation_Summary.md`
- **SLURM Integration**: `docs/SLURM_Aware_Logging_Summary.md`

## 🧪 Testing

**⚠️ Important: Tests should be run after setting up `config/environment.sh` (see [Configuration](#-configuration) section).**

The pipeline includes comprehensive testing:

```bash
# Basic module testing
bash test/test_modules.sh

# Comprehensive functionality testing
bash test/test_comprehensive.sh
```

**Note**: Some tests may show warnings about missing dependencies when functions are tested in isolation. This is expected behavior and does not indicate a problem with the pipeline.

## 📞 Support

For questions or issues:
1. Check the documentation in `docs/`
2. Run the test scripts to verify functionality
3. Review logs for error messages
4. Check SLURM job status and logs

## Acknowledgments

- **Authors**: Phu Khang Ha
- **Contributors**: Paulo Henrique da Silva, Lisa Pasipanodya, Shashi Goonetilleke, Daniel Edge-Garza, Elizabeth Ross
- **Institution**: QAAFI, UQ
- **Pipeline**: Based on GATK Best Practices for variant calling

---

**Version**: 1.0  
**Last Updated**: October 2025  
**License**: Research Use
