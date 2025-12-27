#!/bin/bash
# =============================================================================
# PIPELINE ENVIRONMENT TEMPLATE
# =============================================================================
# Copy this file to config/environment.sh and update values to match your HPC
# environment. The pipeline automatically sources config/environment.sh when
# present; if it is missing, this template is used as a fallback.
# =============================================================================

# -----------------------------------------------------------------------------
# PIPELINE LOCATION
# -----------------------------------------------------------------------------

# Point this to the absolute directory that contains the pipeline checkout.
# $HOME is a sensible default, but the repository can also live on RDM or scratch.
# Always adjust this path to match where you've deployed the code.
PIPELINE_ROOT="${PIPELINE_ROOT:-$HOME/GATK_Pipeline_KH_v1}"

# -----------------------------------------------------------------------------
# STORAGE LOCATIONS
# -----------------------------------------------------------------------------

# Permanent project storage (RDM / QRISdata)
PIPELINE_RDM_BASE="${PIPELINE_RDM_BASE:-/QRISdata/Q8367}"
PIPELINE_RDM_DATASETS="${PIPELINE_RDM_DATASETS:-${PIPELINE_RDM_BASE}/WGS_Reference_Panel}"

# Scratch workspace root for temporary files
PIPELINE_SCRATCH_BASE="${PIPELINE_SCRATCH_BASE:-/scratch/user/uqpha1}"

# Shared project log directory (optional, defaults to scratch if unset)
PIPELINE_LOG_BASE="${PIPELINE_LOG_BASE:-/scratch/user/uqpha1/logs}"

# -----------------------------------------------------------------------------
# REFERENCE DATA
# -----------------------------------------------------------------------------

# Directory containing the reference FASTA and associated index/dict files
PIPELINE_REFERENCE_DIR="/QRISdata/Q8367/Reference_Genome"

# Reference FASTA (including .fai/.dict alongside)
PIPELINE_REFERENCE_FASTA="${PIPELINE_REFERENCE_DIR}/GDDH13_1-1_formatted.fasta"

# Known sites VCF used for BQSR (with .tbi index)
PIPELINE_KNOWN_SITES_VCF="${PIPELINE_RDM_BASE}/Known_Sites/Final_DB_Known-sites_Filtered_by_480K_Unique_HD.vcf"

# Adapter FASTA for Trimmomatic (usually supplied with the module)
PIPELINE_ADAPTER_FASTA="${PIPELINE_RDM_BASE}/Liao_2021_adapter/TruSeq2-PE.fa"

# -----------------------------------------------------------------------------
# SLURM DEFAULTS
# -----------------------------------------------------------------------------

# Default account / project
PIPELINE_SLURM_ACCOUNT="a_qaafi_cas"

# Default partition / queue
PIPELINE_SLURM_PARTITION="general"

# Optional QoS or other flags (leave empty if unused)
PIPELINE_SLURM_QOS=""

# -----------------------------------------------------------------------------
# STEP-SPECIFIC DEFAULTS
# -----------------------------------------------------------------------------

# Step 1A (per-sample calling)
STEP1A_CPUS=10
STEP1A_MEMORY="32G"
STEP1A_TIME="200:00:00"
STEP1A_BATCH_SIZE=50      # Used for optional batching logic

# Step 1B (combine GVCFs)
STEP1B_CPUS=6
STEP1B_MEMORY="36G"
STEP1B_TIME="200:00:00"
STEP1B_ARRAY_LIMIT=25     # Max concurrent array tasks

# Step 1C (Beagle) placeholder defaults
STEP1C_CPUS=8
STEP1C_MEMORY="48G"
STEP1C_TIME="48:00:00"
STEP1C_ACCOUNT=""          # Optional: override default account (leave empty to use PIPELINE_SLURM_ACCOUNT)
STEP1C_PARTITION=""       # Optional: override default partition (leave empty to use PIPELINE_SLURM_PARTITION)
STEP1C_NODES=""           # Optional: override default nodes (leave empty to use PIPELINE_SLURM_NODES)
STEP1C_NTASKS=""          # Optional: override default ntasks (leave empty to use PIPELINE_SLURM_NTASKS)
# Beagle self-imputation toggle: true => Beagle impute=true (self-impute untyped
# markers, no hard-filter); false/empty => phasing-only (impute=false, hard-filter on)
STEP1C_SELF_IMPUTE="false"

# Step 1D (QC analysis)
STEP1D_CPUS=16
STEP1D_MEMORY="256G"
STEP1D_TIME="72:00:00"
STEP1D_ACCOUNT=""          # Optional: override default account (leave empty to use PIPELINE_SLURM_ACCOUNT)
STEP1D_PARTITION=""       # Optional: override default partition (leave empty to use PIPELINE_SLURM_PARTITION)
STEP1D_NODES=""           # Optional: override default nodes (leave empty to use PIPELINE_SLURM_NODES)
STEP1D_NTASKS=""          # Optional: override default ntasks (leave empty to use PIPELINE_SLURM_NTASKS)
# Step 1D behavior toggles
STEP1D_PCA_ONLY="false"            # true => run PCA only
STEP1D_REMOVE_RELATIVES="false"    # true => drop relatives before PCA/QC
STEP1D_PCA_SHOW_LABELS="true"      # show sample labels on PCA plots
STEP1D_PCA_USE_GGREPEL="true"      # use ggrepel to avoid label overlap
STEP1D_PCA_LABEL_SIZE="1.5"        # PCA label size
STEP1D_PCA_DIR="pca_analysis"      # output directory name for PCA assets
STEP1D_PCA_MERGED_PATTERN="*merged*.vcf.gz,*merge*.vcf.gz"  # comma-separated glob(s) to locate merged VCFs for PCA
STEP1D_PCA_FORCE_CONCAT="false"    # true => always concatenate per-chrom VCFs (ignore merged detection)
STEP1D_PCA_MERGED_EXCLUDE_CHR="true"  # true => ignore merged candidates containing "Chr" in filename
STEP1D_DUPLICATE_MODE="flag"       # off|flag|remove; KING-based duplicate detection for PCA
STEP1D_DUPLICATE_KING_THRESHOLD="0.45"  # kinship threshold for duplicate detection (e.g. 0.45)
STEP1D_AF_PLOTS_DIR="af_distribution_plots"  # output directory name for allele frequency plots
STEP1D_AF_HIST_BINS="50"            # histogram bins for allele frequency plots

# -----------------------------------------------------------------------------
# ADVANCED / OPTIONAL SETTINGS
# -----------------------------------------------------------------------------

# Logging and global behavior
PIPELINE_LOG_LEVEL="INFO"
LOG_TO_CONSOLE="true"
LOG_TO_FILE="true"
PIPELINE_VERBOSE_LOGGING="false"
PIPELINE_GLOBAL_DRY_RUN="false"

# Monitoring / notifications
MONITOR_INTERVAL="30"        # seconds between monitor polls
MAX_MONITOR_ATTEMPTS="100"   # max monitor iterations
EMAIL_NOTIFICATIONS="false"  # set true to send job notifications
EMAIL_ADDRESS=""             # email target when notifications enabled

# Resume / cleanup / backups
ENABLE_RESUME="true"
CLEANUP_TEMP_FILES="true"
CLEANUP_INTERMEDIATE_FILES="false"
PIPELINE_ENABLE_BACKUP="true"
PIPELINE_BACKUP_STEPS="3,4,5,6"
PIPELINE_SHARED_REF_TIMEOUT="120"
PIPELINE_SHARED_REF_POLL_INTERVAL="3"

# Runtime / tool knobs
GATK_BATCH_SIZE="50"
GATK_READER_THREADS="4"
BWA_MEM_FLAGS="-M"
