#!/bin/bash
################################################################################
# Master VCF Analysis Pipeline
# =============================================================================
# This script orchestrates the complete VCF analysis workflow:
#   1. Extract mean depth TSV from VCF files (if not exists)
#   2. Generate per-chromosome depth plots
#   3. Generate per-site missingness vs depth plots
#   4. Combine all missingness plots into single image
#
# Usage:
#   Interactive: bash master_vcf_analysis.sh
#   SLURM: sbatch master_vcf_analysis.sh
#
#
# NOTE: This script always filters inputs to biallelic SNPs (bcftools view -m2 -M2 -v snps)
# to ensure consistent QC metrics across datasets.
################################################################################

set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODULE_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
if [ -n "${PIPELINE_ROOT:-}" ]; then
    if [ -d "${PIPELINE_ROOT}" ]; then
        PIPELINE_ROOT="$(cd "${PIPELINE_ROOT}" && pwd)"
    else
        echo "[master_vcf_analysis] ⚠️  Provided PIPELINE_ROOT (${PIPELINE_ROOT}) does not exist; falling back to template-relative path." >&2
        PIPELINE_ROOT="$(cd "${MODULE_DIR}/.." && pwd)"
    fi
else
    PIPELINE_ROOT="$(cd "${MODULE_DIR}/.." && pwd)"
fi
export PIPELINE_ROOT
if [ -f "${PIPELINE_ROOT}/config/pipeline_config.sh" ]; then
    source "${PIPELINE_ROOT}/config/pipeline_config.sh"
fi
if [ -f "${PIPELINE_ROOT}/lib/logging.sh" ]; then
    source "${PIPELINE_ROOT}/lib/logging.sh"
    init_logging "step1d" "pipeline" "${DATASET_NAME:-${DATASET:-}}"
fi
trap 'log_error "Unexpected failure in master_vcf_analysis.sh"; exit 1' ERR

usage() {
    cat <<'EOF'
Usage: master_vcf_analysis.sh [--dry-run] [--beagle] [--qc] [--PCA] [--duplicate-check] [--remove-relatives] [--help]

Options:
  --dry-run, -n      Print the actions that would be performed without creating files
  --beagle           Treat input VCFs as Beagle-imputed (expected INFO tags: AF, DR2, IMP)
  --qc               Run QC-only mode (metrics + plots); default if no mode is specified
  --PCA              Run PCA-only mode (expects merged or per-chrom VCFs)
  --duplicate-check  Run KING-based duplicate detection only (writes duplicate pairs/list)
  --remove-relatives Remove close relatives before PCA (requires --PCA)
  --help, -h         Show this message and exit
EOF
}

DRY_RUN=false
BEAGLE_MODE=false
REMOVE_RELATIVES_CLI=false
SHOW_USAGE=false
EXIT_ERROR=false
MODE_QC=false
MODE_PCA=false
MODE_DUP_CHECK=false

for arg in "$@"; do
    case "$arg" in
        --dry-run|-n)
            DRY_RUN=true
            ;;
        --beagle)
            BEAGLE_MODE=true
            ;;
        --qc)
            MODE_QC=true
            ;;
        --PCA|--pca)
            MODE_PCA=true
            ;;
        --duplicate-check)
            MODE_DUP_CHECK=true
            ;;
        --remove-relatives)
            REMOVE_RELATIVES_CLI=true
            ;;
        --help|-h)
            SHOW_USAGE=true
            ;;
        *)
            echo "❌ Unknown option: ${arg}" >&2
            SHOW_USAGE=true
            EXIT_ERROR=true
            ;;
    esac
done

if [ "${SHOW_USAGE}" = "true" ]; then
    usage
    if [ "${EXIT_ERROR}" = "true" ]; then
        exit 1
    else
        exit 0
    fi
fi

mode_count=0
${MODE_QC} && mode_count=$((mode_count + 1))
${MODE_PCA} && mode_count=$((mode_count + 1))
${MODE_DUP_CHECK} && mode_count=$((mode_count + 1))

if [ "${mode_count}" -gt 1 ]; then
    echo "❌ Please select only one mode: --qc, --PCA, or --duplicate-check." >&2
    exit 1
fi

if [ "${mode_count}" -eq 0 ]; then
    MODE_QC=true
fi

if [ "${REMOVE_RELATIVES_CLI}" = "true" ] && [ "${MODE_PCA}" != "true" ]; then
    echo "❌ --remove-relatives requires --PCA." >&2
    exit 1
fi

# Utility helpers used during configuration
normalize_bool() {
    local value="$1"
    if [ -z "${value}" ]; then
        echo "false"
        return
    fi
    value="$(printf '%s' "${value}" | tr '[:upper:]' '[:lower:]')"
    case "${value}" in
        true|1|yes|y)
            echo "true"
            ;;
        *)
            echo "false"
            ;;
    esac
}

# =============================================================================
# CONFIGURATION (simple, HPC-friendly)
# - You can override any of these via environment variables when calling sbatch
#   e.g. sbatch --export=ALL,VCF_DIR=/scratch/me/vcfs master_vcf_analysis.sh
# =============================================================================

# VCF file directory (default: current directory)
VCF_DIR="${VCF_DIR:-$PWD}"

# VCF file naming pattern (00-17 → Chr00.vcf.gz, Chr01.vcf.gz, ...)
VCF_PATTERN="${VCF_PATTERN:-Chr%02d.vcf.gz}"

# Working directory (default: VCF_DIR)
WORK_DIR="${WORK_DIR:-${VCF_DIR}}"

# R scripts directory (default: modules/step1d/Rscripts)
R_SCRIPTS_DIR="${R_SCRIPTS_DIR:-${MODULE_DIR}/Rscripts}"

# Output datasets and directories
MEAN_DP_TSV="${MEAN_DP_TSV:-SNP_site_meanDP.tsv}"
SITE_METRICS_TSV="${SITE_METRICS_TSV:-variant_site_metrics.tsv}"
METRICS_BY_CHROM_DIR="${METRICS_BY_CHROM_DIR:-site_metrics_per_chromosome}"
DEPTH_PLOTS_DIR="${DEPTH_PLOTS_DIR:-depth_plots}"
MISSINGNESS_PLOTS_DIR="${MISSINGNESS_PLOTS_DIR:-missingness_plots}"
DEPTH_MISSINGNESS_DIR="${DEPTH_MISSINGNESS_DIR:-depth_vs_missingness}"
SITE_QUALITY_DIR="${SITE_QUALITY_DIR:-site_quality_plots}"
HETEROZYGOSITY_DIR="${HETEROZYGOSITY_DIR:-heterozygosity_plots}"
QD_PLOTS_DIR="${QD_PLOTS_DIR:-quality_by_depth_plots}"
CALL_RATE_HEATMAP_DIR="${CALL_RATE_HEATMAP_DIR:-call_rate_heatmaps}"
PLOT_IMAGE_FORMAT="${PLOT_IMAGE_FORMAT:-png}"
CALL_RATE_HEATMAP_BINS="${CALL_RATE_HEATMAP_BINS:-100}"
MEAN_DP_PATH="${WORK_DIR}/${MEAN_DP_TSV}"
SITE_METRICS_PATH="${WORK_DIR}/${SITE_METRICS_TSV}"
METRICS_BY_CHROM_PATH="${WORK_DIR}/${METRICS_BY_CHROM_DIR}"
DEPTH_OUTPUT_DIR="${WORK_DIR}/${DEPTH_PLOTS_DIR}"
MISSINGNESS_OUTPUT_DIR="${WORK_DIR}/${MISSINGNESS_PLOTS_DIR}"
DEPTH_MISS_OUTPUT_DIR="${WORK_DIR}/${DEPTH_MISSINGNESS_DIR}"
SITE_QUALITY_OUTPUT_DIR="${WORK_DIR}/${SITE_QUALITY_DIR}"
HETEROZYGOSITY_OUTPUT_DIR="${WORK_DIR}/${HETEROZYGOSITY_DIR}"
QD_OUTPUT_DIR="${WORK_DIR}/${QD_PLOTS_DIR}"
QD_OUTPUT_FILE="${QD_OUTPUT_DIR}/quality_by_depth_hist.${PLOT_IMAGE_FORMAT}"
CALL_RATE_OUTPUT_DIR="${WORK_DIR}/${CALL_RATE_HEATMAP_DIR}"
AF_PLOTS_DIR="${AF_PLOTS_DIR:-${STEP1D_AF_PLOTS_DIR:-af_distribution_plots}}"
AF_HIST_BINS="${AF_HIST_BINS:-${STEP1D_AF_HIST_BINS:-50}}"
AF_OUTPUT_DIR="${WORK_DIR}/${AF_PLOTS_DIR}"
HAS_DP=false
HAS_QD=false
PCA_DIR_NAME="${STEP1D_PCA_DIR:-pca_analysis}"
PCA_OUTPUT_DIR="${WORK_DIR}/${PCA_DIR_NAME}"
RUN_QC="${MODE_QC}"
RUN_PCA="${MODE_PCA}"
RUN_DUP_CHECK="${MODE_DUP_CHECK}"
REMOVE_RELATIVES="${STEP1D_REMOVE_RELATIVES:-false}"
PLINK2_BIN_PATH="${PLINK2_BIN:-plink2}"
BCFTOOLS_BIN_PATH="${BCFTOOLS_BIN:-bcftools}"
PCA_MERGED_PATTERN="${STEP1D_PCA_MERGED_PATTERN:-*merged*.vcf.gz,*merge*.vcf.gz}"
PCA_FORCE_CONCAT="${STEP1D_PCA_FORCE_CONCAT:-false}"
PCA_MERGED_EXCLUDE_CHR="${STEP1D_PCA_MERGED_EXCLUDE_CHR:-true}"
DUPLICATE_MODE="${STEP1D_DUPLICATE_MODE:-flag}"
DUPLICATE_KING_THRESHOLD="${STEP1D_DUPLICATE_KING_THRESHOLD:-0.45}"
SHOW_LABELS="${STEP1D_PCA_SHOW_LABELS:-true}"
LABEL_SIZE="${STEP1D_PCA_LABEL_SIZE:-1.5}"
USE_GGREPEL="${STEP1D_PCA_USE_GGREPEL:-true}"
export STEP1D_PCA_FORCE_CONCAT="${PCA_FORCE_CONCAT}"
export STEP1D_PCA_MERGED_EXCLUDE_CHR="${PCA_MERGED_EXCLUDE_CHR}"
PCA_EXECUTED="false"
PCA_SKIP_REASON=""

REMOVE_RELATIVES="$(normalize_bool "${REMOVE_RELATIVES}")"
if [ "${REMOVE_RELATIVES_CLI}" = "true" ]; then
    REMOVE_RELATIVES="true"
fi

if [ "${REMOVE_RELATIVES}" = "true" ] && [ "${RUN_PCA}" != "true" ]; then
    log_error "--remove-relatives requires --PCA."
    exit 1
fi

declare -a VCF_TARGETS=()
if [ -n "${VCF_INCLUDE_FILENAMES:-}" ]; then
    read -r -a VCF_TARGETS <<< "${VCF_INCLUDE_FILENAMES}"
else
    for chr_num in $(seq 0 17); do
        VCF_TARGETS+=("$(printf "${VCF_PATTERN}" "${chr_num}")")
    done
fi

# Limit targets to VCFs that actually exist in VCF_DIR
declare -a EXISTING_VCF_TARGETS=()
for vcf_basename in "${VCF_TARGETS[@]}"; do
    if [ -f "${VCF_DIR}/${vcf_basename}" ]; then
        EXISTING_VCF_TARGETS+=("${vcf_basename}")
    fi
done
VCF_TARGETS=("${EXISTING_VCF_TARGETS[@]}")

if [ ${#VCF_TARGETS[@]} -eq 0 ]; then
    log_error "No VCF files specified for analysis."
    exit 1
fi

declare -a CHR_LABELS=()
declare -A CHR_TO_VCF=()
for vcf_basename in "${VCF_TARGETS[@]}"; do
    chr_name="${vcf_basename%.vcf.gz}"
    chr_name="${chr_name%.vcf}"
    CHR_LABELS+=("${chr_name}")
    CHR_TO_VCF["${chr_name}"]="${vcf_basename}"
done

EXPECTED_CHROM_COUNT=${#CHR_LABELS[@]}

# Conda environment name for R (optional; if not present, will use PATH)
CONDA_ENV="${CONDA_ENV:-rplot}"

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

log_section() {
    echo ""
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "$1"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo ""
}

log_success() { log_info "$1"; }
log_warning() { log_warn "$1"; }

check_file() {
    if [ ! -f "$1" ]; then
        log_error "Required file not found: $1"
        return 1
    fi
    return 0
}

check_command() {
    if ! command -v "$1" &> /dev/null; then
        log_error "Required command not found: $1"
        return 1
    fi
    return 0
}

ensure_r_packages() {
    local missing=()
    for pkg in "$@"; do
        if ! Rscript -e "quit(status = ifelse(requireNamespace('${pkg}', quietly = TRUE), 0, 1))" >/dev/null 2>&1; then
            missing+=("${pkg}")
        fi
    done
    if [ ${#missing[@]} -gt 0 ]; then
        log_error "Missing required R packages: ${missing[*]}"
        log_error "Install them inside the '${CONDA_ENV}' environment, e.g.:"
        log_error ""
        for pkg in "${missing[@]}"; do
            log_error "    conda install -n ${CONDA_ENV} r-${pkg}"
        done
        log_error ""
        log_error "Or activate ${CONDA_ENV} and run:"
        log_error "    Rscript -e \"install.packages(c('${missing[*]}'), repos='https://cloud.r-project.org')\""
        exit 1
    fi
}

# =============================================================================
# SETUP
# =============================================================================

log_section "Master VCF Analysis Pipeline - Starting"

log_info "Pipeline started"
log_info "VCF Directory: ${VCF_DIR}"
log_info "Working Directory: ${WORK_DIR}"
log_info "R Scripts Directory: ${R_SCRIPTS_DIR}"
if [ "${BEAGLE_MODE}" = "true" ]; then
    log_info "Beagle mode enabled: using AF/DR2/IMP metrics and skipping depth-dependent outputs."
fi
log_info "Selected mode -> QC: ${RUN_QC}, PCA: ${RUN_PCA}, Duplicate-check: ${RUN_DUP_CHECK}"
if [ "${RUN_PCA}" = "true" ]; then
    log_info "PCA settings: remove relatives: ${REMOVE_RELATIVES}; duplicate mode: ${DUPLICATE_MODE}; KING threshold: ${DUPLICATE_KING_THRESHOLD}"
fi

if [ "${DRY_RUN}" = "true" ]; then
    log_info "Dry-run mode enabled. No files will be created."
fi

# Change to working directory
cd "${WORK_DIR}" || {
    log_error "Cannot change to working directory: ${WORK_DIR}"
    exit 1
}

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

BCFTOOLS_MODULE="${BCFTOOLS_MODULE:-bcftools/1.18-GCC-12.3.0}"
if module load "${BCFTOOLS_MODULE}" >/dev/null 2>&1; then
    log_success "Loaded ${BCFTOOLS_MODULE} module"
else
    log_error "Failed to load ${BCFTOOLS_MODULE} module"
    exit 1
fi

if [ "${NEED_PLINK}" = "true" ]; then
    if module load plink/2.00a3.6-gcc-11.3.0 >/dev/null 2>&1; then
        log_success "Loaded plink/2.00a3.6-gcc-11.3.0 module"
    else
        log_error "Failed to load plink/2.00a3.6-gcc-11.3.0 module"
        exit 1
    fi
else
    log_info "Skipping plink module load (not required for this mode)"
fi

if [ "${NEED_R}" = "true" ]; then
    if [ -f "$ROOTMINIFORGE/etc/profile.d/conda.sh" ]; then
        source "$ROOTMINIFORGE/etc/profile.d/conda.sh"
        if conda activate "${CONDA_ENV}" >/dev/null 2>&1; then
            log_success "Conda environment '${CONDA_ENV}' activated"
        else
            log_error "Failed to activate conda environment: ${CONDA_ENV}"
            exit 1
        fi
    else
        log_warning "Conda initialization script not found; assuming R is available in PATH"
    fi
fi

# Check required commands
log_info "Checking required commands..."
check_command bcftools || exit 1
if [ "${NEED_PLINK}" = "true" ]; then
    check_command plink2 || exit 1
fi
if [ "${NEED_R}" = "true" ]; then
    check_command Rscript || exit 1
    ensure_r_packages ggplot2 data.table ragg scales
fi
log_success "Required commands and packages available"

# =============================================================================
# STEP 1: BUILD SITE-LEVEL METRICS DATASET
# =============================================================================

if [ "${RUN_QC}" = "true" ]; then

log_section "Step 1: Build Site-Level Metrics Dataset"
if [ "${BEAGLE_MODE}" = "true" ]; then
    SITE_METRICS_HEADER=$'CHROM\tPOS\tQUAL\tAF\tDR2\tIMP\tCALL_RATE\tMISSING_RATE\tHETEROZYGOUS_RATE\tCALLED_GENOTYPES\tMISSING_GENOTYPES\tTOTAL_GENOTYPES\tHETEROZYGOUS_COUNT'
    HAS_DP=false
    HAS_QD=false
else
    SITE_METRICS_HEADER=$'CHROM\tPOS\tQUAL\tQD\tAC\tAF\tINBREEDING_COEFF\tEXCESS_HET\tMQ\tMEAN_DEPTH\tCALL_RATE\tMISSING_RATE\tHETEROZYGOUS_RATE\tDP_NON_MISSING\tCALLED_GENOTYPES\tMISSING_GENOTYPES\tTOTAL_GENOTYPES\tHETEROZYGOUS_COUNT'
    HAS_DP=true
    HAS_QD=true
fi

if [ "${DRY_RUN}" = "true" ]; then
    log_info "[dry-run] Would ensure directory exists: ${METRICS_BY_CHROM_PATH}"
else
    mkdir -p "${METRICS_BY_CHROM_PATH}"
fi

REGENERATE_METRICS=false
if [ -f "${SITE_METRICS_PATH}" ]; then
    header_line="$(head -n 1 "${SITE_METRICS_PATH}" 2>/dev/null || true)"
    if [ -z "${header_line}" ] || ! printf '%s' "${header_line}" | grep -q $'\t'; then
        log_warning "Existing site metrics TSV appears malformed (missing tab-separated header). Regenerating metrics."
        REGENERATE_METRICS=true
    else
        if [ "${BEAGLE_MODE}" = "true" ]; then
            for marker in AF DR2 IMP; do
                if ! printf '%s\n' "${header_line}" | tr '\t' '\n' | grep -Fxq "${marker}"; then
                    log_warning "Existing site metrics missing '${marker}' column required for --beagle mode. Regenerating metrics."
                    REGENERATE_METRICS=true
                    break
                fi
            done
            if [ "${REGENERATE_METRICS}" = false ] && printf '%s\n' "${header_line}" | tr '\t' '\n' | grep -Fxq "MEAN_DEPTH"; then
                log_warning "Existing site metrics include depth columns but --beagle mode requested. Regenerating metrics."
                REGENERATE_METRICS=true
            fi
        else
            for marker in MEAN_DEPTH QD MQ; do
                if ! printf '%s\n' "${header_line}" | tr '\t' '\n' | grep -Fxq "${marker}"; then
                    log_warning "Existing site metrics missing '${marker}' column. Regenerating metrics."
                    REGENERATE_METRICS=true
                    break
                fi
            done
            if [ "${REGENERATE_METRICS}" = false ] && printf '%s\n' "${header_line}" | tr '\t' '\n' | grep -Fxq "DR2"; then
                log_warning "Existing site metrics appear to have been generated in --beagle mode. Regenerating metrics for full dataset."
                REGENERATE_METRICS=true
            fi
        fi
        if [ "${REGENERATE_METRICS}" = false ]; then
            existing_metrics=$(find "${METRICS_BY_CHROM_PATH}" -maxdepth 1 -type f -name '*_metrics.tsv' | awk 'END{print NR}')
            if [ "${existing_metrics}" -ne "${EXPECTED_CHROM_COUNT}" ]; then
                log_warning "Per-chromosome metric files (${existing_metrics}) do not match expected count (${EXPECTED_CHROM_COUNT}). Skipping regeneration since combined TSV looks valid."
            fi
        fi
    fi
fi

if [ "${REGENERATE_METRICS}" = true ]; then
    rm -f "${SITE_METRICS_PATH}"
    find "${METRICS_BY_CHROM_PATH}" -maxdepth 1 -type f -name '*_metrics.tsv' -delete 2>/dev/null || true
    find "${WORK_DIR}" -maxdepth 1 -type f -name '*_snps.vcf.gz*' -delete 2>/dev/null || true
fi

if [ -f "${SITE_METRICS_PATH}" ]; then
    log_success "Site metrics TSV already exists: ${SITE_METRICS_PATH}"
    log_info "Skipping regeneration of site metrics"
else
    if [ "${DRY_RUN}" = "true" ]; then
        log_info "[dry-run] Target VCF files: ${VCF_TARGETS[*]}"
        log_info "[dry-run] Would create filtered SNP VCFs at ${WORK_DIR}/<chromosome>_snps.vcf.gz"
        log_info "[dry-run] Would generate per-chromosome metrics in ${METRICS_BY_CHROM_PATH}/<chromosome>_metrics.tsv"
        log_info "[dry-run] Would combine metrics into ${SITE_METRICS_PATH}"
    else
        log_info "Generating site metrics TSV from per-chromosome VCFs..."
        
        if [ -z "${TMPDIR:-}" ]; then
            log_error "TMPDIR is not defined. Please ensure the job allocates scratch space."
            exit 1
        fi
        TEMP_DIR=$(mktemp -d "${TMPDIR%/}/step1d_site_metrics.XXXXXX")
        log_info "Created temporary directory: ${TEMP_DIR}"
        
        declare -a METRIC_FILES=()
        declare -a NEW_METRIC_FILES=()
        
        for vcf_basename in "${VCF_TARGETS[@]}"; do
            chr_name="${vcf_basename%.vcf.gz}"
            chr_name="${chr_name%.vcf}"
            vcf_file="${VCF_DIR}/${vcf_basename}"
            filtered_temp="${TEMP_DIR}/${chr_name}_snps.vcf.gz"
            filtered_output="${WORK_DIR}/${chr_name}_snps.vcf.gz"
            metrics_output="${METRICS_BY_CHROM_PATH}/${chr_name}_metrics.tsv"
            tmp_metrics="${TEMP_DIR}/${chr_name}_metrics.tsv"
            
            if [ ! -f "${vcf_file}" ]; then
                log_warning "${chr_name}: VCF not found (${vcf_file}), skipping"
                continue
            fi

            if [ -f "${metrics_output}" ] && [ "${REGENERATE_METRICS}" = false ]; then
                log_info "${chr_name}: metrics TSV already present (${metrics_output}); skipping regeneration."
                if [ ! -f "${filtered_output}" ]; then
                    log_warning "${chr_name}: filtered SNP VCF not found at ${filtered_output} (metrics will be reused)."
                fi
                METRIC_FILES+=("${metrics_output}")
                continue
            fi
            
            log_info "${chr_name}: filtering to biallelic SNPs"
            if ! bcftools view -m2 -M2 -v snps -Oz -o "${filtered_temp}" "${vcf_file}"; then
                log_error "${chr_name}: failed to filter biallelic SNPs"
                rm -rf "${TEMP_DIR}"
                exit 1
            fi
            if ! bcftools index -t "${filtered_temp}"; then
                log_error "${chr_name}: failed to index filtered SNP VCF"
                rm -rf "${TEMP_DIR}"
                exit 1
            fi
            mv -f "${filtered_temp}" "${filtered_output}"
            mv -f "${filtered_temp}.tbi" "${filtered_output}.tbi"
            log_info "${chr_name}: filtered SNP VCF written to ${filtered_output}"

            log_info "${chr_name}: extracting site metrics..."

            if [ "${BEAGLE_MODE}" = "true" ]; then
                query_format='%CHROM\t%POS\t%QUAL\t%INFO/AF\t%INFO/DR2\t%INFO/IMP\t[%GT\t]\n'
                if ! bcftools query -f "${query_format}" "${filtered_output}" | \
                    awk -v OFS='\t' -v header="${SITE_METRICS_HEADER}" '
                        BEGIN {
                            print header;
                        }
                        {
                            sub(/\t$/, "", $0);
                            chrom = $1;
                            pos = $2;
                            qual = $3;
                            af = $4;
                            dr2 = $5;
                            imp = $6;
                            total_genotypes = NF - 6;
                            called = 0;
                            missing = 0;
                            hetero = 0;

                            for (i = 7; i <= NF; i++) {
                                gt = $i;
                                if (gt == "" || gt == "." || gt == "./." || gt == ".|.") {
                                    missing++;
                                } else {
                                    called++;
                                    split(gt, alleles, /[\/|]/);
                                    if (length(alleles) >= 2) {
                                        a1 = alleles[1];
                                        a2 = alleles[2];
                                        if (a1 != "." && a2 != "." && a1 != a2) {
                                            hetero++;
                                        }
                                    }
                                }
                            }

                            qual_value = (qual == "." ? "NA" : qual);
                            af_value = (af == "" || af == "." ? "NA" : af);
                            dr2_value = (dr2 == "" || dr2 == "." ? "NA" : dr2);
                            imp_value = (imp == "" || imp == "." ? "0" : "1");

                            total_genotypes = (total_genotypes > 0 ? total_genotypes : 0);
                            call_rate = (total_genotypes > 0 ? called / total_genotypes : 0);
                            missing_rate = (total_genotypes > 0 ? missing / total_genotypes : 0);
                            hetero_rate = (called > 0 ? hetero / called : 0);

                            printf "%s\t%s\t%s\t%s\t%s\t%s\t%.6f\t%.6f\t%.6f\t%d\t%d\t%d\t%d\n",
                                   chrom, pos, qual_value, af_value, dr2_value, imp_value,
                                   call_rate, missing_rate, hetero_rate, called, missing, total_genotypes, hetero;
                        }
                    ' > "${tmp_metrics}"
                then
                    log_error "${chr_name}: failed to generate site metrics"
                    rm -rf "${TEMP_DIR}"
                    exit 1
                fi
            else
                if ! bcftools query \
                    -f '%CHROM\t%POS\t%QUAL\t%INFO/QD\t%INFO/AC\t%INFO/AF\t%INFO/InbreedingCoeff\t%INFO/ExcessHet\t%INFO/MQ\t[%GT:%DP\t]\n' \
                    "${filtered_output}" | \
                    awk -v OFS='\t' -v header="${SITE_METRICS_HEADER}" '
                        BEGIN {
                            print header;
                        }
                        {
                            sub(/\t$/, "", $0);
                            chrom = $1;
                            pos = $2;
                            qual = $3;
                            qd = $4;
                            ac = $5;
                            af = $6;
                            inbreeding = $7;
                            excesshet = $8;
                            mq = $9;
                            sum_dp = 0;
                            dp_non_missing = 0;
                            total = 0;
                            missing = 0;
                            hetero = 0;
                            called = 0;

                            qd_value = (qd == "." || qd == "" || qd == "NA") ? "NA" : qd;
                            ac_value = (ac == "" || ac == "NA" || ac == ".") ? "NA" : ac;
                            af_value = (af == "" || af == "NA" || af == ".") ? "NA" : af;
                            inbreed_value = (inbreeding == "" || inbreeding == "NA" || inbreeding == ".") ? "NA" : inbreeding;
                            excesshet_value = (excesshet == "" || excesshet == "NA" || excesshet == ".") ? "NA" : excesshet;
                            mq_value = (mq == "" || mq == "NA" || mq == ".") ? "NA" : mq;

                            for (i = 10; i <= NF; i++) {
                                field = $i;
                                if (field == "") {
                                    continue;
                                }

                                total++;
                                split(field, parts, ":");
                                gt = parts[1];
                                dp = (length(parts) > 1 ? parts[2] : "");

                                if (dp != "" && dp != "." && dp != "./." && dp != ".|." && dp ~ /^-?[0-9]+(\.[0-9]+)?$/) {
                                    sum_dp += dp;
                                    dp_non_missing++;
                                }

                                if (gt != "" && gt != "." && gt != "./." && gt != ".|.") {
                                    split(gt, alleles, /[\/|]/);
                                    if (length(alleles) >= 2) {
                                        a1 = alleles[1];
                                        a2 = alleles[2];
                                        if (a1 != "." && a2 != ".") {
                                            called++;
                                            if (a1 != a2) {
                                                hetero++;
                                            }
                                        }
                                    }
                                } else {
                                    missing++;
                                }
                            }

                            mean_depth = (dp_non_missing > 0 ? sum_dp / dp_non_missing : 0);
                            call_rate = (total > 0 ? (total - missing) / total : 0);
                            missing_rate = 1 - call_rate;
                            hetero_rate = (called > 0 ? hetero / called : 0);
                            qual_value = (qual == "." ? "NA" : qual);

                            printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%.6f\t%.6f\t%.6f\t%.6f\t%d\t%d\t%d\t%d\t%d\n",
                                   chrom, pos, qual_value, qd_value, ac_value, af_value, inbreed_value, excesshet_value, mq_value,
                                   mean_depth, call_rate, missing_rate, hetero_rate, dp_non_missing, called, missing, total, hetero;
                        }
                    ' > "${tmp_metrics}"
                then
                    log_error "${chr_name}: failed to generate site metrics"
                    rm -rf "${TEMP_DIR}"
                    exit 1
                fi
            fi

            if [ ! -s "${tmp_metrics}" ]; then
                log_warning "${chr_name}: no variant records found, skipping"
                rm -f "${tmp_metrics}"
                continue
            fi

            mv "${tmp_metrics}" "${metrics_output}"
            log_success "${chr_name}: metrics written to ${metrics_output}"
            METRIC_FILES+=("${metrics_output}")
            NEW_METRIC_FILES+=("${metrics_output}")
        done
        
        if [ ${#METRIC_FILES[@]} -eq 0 ]; then
            log_error "No per-chromosome metrics files are available to combine. Aborting."
            rm -rf "${TEMP_DIR}"
            exit 1
        fi

        if [ ${#NEW_METRIC_FILES[@]} -gt 0 ]; then
            log_success "Generated ${#NEW_METRIC_FILES[@]} per-chromosome metrics file(s)."
        fi
        reused_count=$(( ${#METRIC_FILES[@]} - ${#NEW_METRIC_FILES[@]} ))
        if [ ${reused_count} -gt 0 ]; then
            log_info "Reusing ${reused_count} existing per-chromosome metrics file(s)."
        fi

        SITE_METRICS_TEMP="${TEMP_DIR}/combined_metrics.tsv"
        printf '%s\n' "${SITE_METRICS_HEADER}" > "${SITE_METRICS_TEMP}"
        for metrics_file in "${METRIC_FILES[@]}"; do
            if [ ! -f "${metrics_file}" ]; then
                log_error "Expected metrics file missing: ${metrics_file}"
                rm -rf "${TEMP_DIR}"
                exit 1
            fi
            tail -n +2 "${metrics_file}" >> "${SITE_METRICS_TEMP}"
        done
        
        mv "${SITE_METRICS_TEMP}" "${SITE_METRICS_PATH}"
        rm -rf "${TEMP_DIR}"
        
        log_success "Site metrics TSV created: ${SITE_METRICS_PATH}"
        DATA_LINES=$(($(wc -l < "${SITE_METRICS_PATH}") - 1))
        log_info "Variant sites recorded: ${DATA_LINES}"
    fi
fi

if [ "${HAS_DP}" = "true" ]; then
    if [ "${DRY_RUN}" = "true" ]; then
        if [ -f "${MEAN_DP_PATH}" ]; then
            LINES=$(wc -l < "${MEAN_DP_PATH}")
            SIZE=$(du -h "${MEAN_DP_PATH}" | cut -f1)
            log_success "Mean depth TSV already present: ${MEAN_DP_PATH}"
            log_info "File size: ${SIZE}, Lines: ${LINES}"
        else
            log_info "[dry-run] Would derive mean depth TSV from ${SITE_METRICS_PATH} → ${MEAN_DP_PATH}"
        fi
    else
        if [ -f "${MEAN_DP_PATH}" ]; then
            LINES=$(wc -l < "${MEAN_DP_PATH}")
            SIZE=$(du -h "${MEAN_DP_PATH}" | cut -f1)
            log_success "Mean depth TSV available: ${MEAN_DP_PATH}"
            log_info "File size: ${SIZE}, Lines: ${LINES}"
        else
            log_info "Deriving mean depth TSV from site metrics..."
            {
                printf "CHROM\tPOS\tMEAN_DEPTH\n"
                awk 'NR == 1 { next } { printf "%s\t%s\t%s\n", $1, $2, $10 }' "${SITE_METRICS_PATH}"
            } > "${MEAN_DP_PATH}"
            
            if [ -f "${MEAN_DP_PATH}" ]; then
                LINES=$(wc -l < "${MEAN_DP_PATH}")
                SIZE=$(du -h "${MEAN_DP_PATH}" | cut -f1)
                log_success "Mean depth TSV created: ${MEAN_DP_PATH}"
                log_info "File size: ${SIZE}, Lines: ${LINES}"
            else
                log_error "Failed to derive mean depth TSV from site metrics"
                exit 1
            fi
        fi
    fi
else
    if [ -f "${MEAN_DP_PATH}" ]; then
        log_warning "Skipping mean depth TSV; existing file left untouched (metric unavailable in current schema)."
    else
        log_info "Skipping mean depth TSV (metric unavailable in current schema)."
    fi
fi

# =============================================================================
# STEP 2: ALLELE FREQUENCY DISTRIBUTION PLOTS
# =============================================================================

log_section "Step 2: Allele Frequency Distribution Plots"

AF_R_SCRIPT="${R_SCRIPTS_DIR}/plot_af_distribution.R"

if ! check_file "${AF_R_SCRIPT}"; then
    log_error "R script not found: ${AF_R_SCRIPT}"
    exit 1
fi

if [ "${DRY_RUN}" = "true" ]; then
    log_info "[dry-run] Would create directory: ${AF_OUTPUT_DIR}"
    log_info "[dry-run] Would run: Rscript ${AF_R_SCRIPT} ${SITE_METRICS_PATH} ${AF_OUTPUT_DIR} ${PLOT_IMAGE_FORMAT} ${AF_HIST_BINS}"
else
    mkdir -p "${AF_OUTPUT_DIR}"
    existing_af_plots=$(find "${AF_OUTPUT_DIR}" -maxdepth 1 -type f -name "*_af_distribution.${PLOT_IMAGE_FORMAT}" | awk 'END{print NR}')
    if [ "${existing_af_plots}" -gt 0 ]; then
        log_info "Allele frequency plots already exist; regenerating to ensure consistency."
    fi
    if Rscript "${AF_R_SCRIPT}" "${SITE_METRICS_PATH}" "${AF_OUTPUT_DIR}" "${PLOT_IMAGE_FORMAT}" "${AF_HIST_BINS}"; then
        AF_COUNT=$(find "${AF_OUTPUT_DIR}" -maxdepth 1 -type f -name "*_af_distribution.${PLOT_IMAGE_FORMAT}" | awk 'END{print NR}')
        log_success "Generated ${AF_COUNT} allele frequency plots"
        log_info "Location: ${AF_OUTPUT_DIR}/"
    else
        log_error "Allele frequency plotting script failed"
        exit 1
    fi
fi

# =============================================================================
# STEP 3: GENERATE MEAN DEPTH PLOTS
# =============================================================================

log_section "Step 3: Generate Mean Depth Plots"

if [ "${HAS_DP}" = "true" ]; then
    DEPTH_R_SCRIPT="${R_SCRIPTS_DIR}/plot_depth_vs_position.R"

    if ! check_file "${DEPTH_R_SCRIPT}"; then
        log_error "R script not found: ${DEPTH_R_SCRIPT}"
        exit 1
    fi

    if [ "${DRY_RUN}" = "true" ]; then
        log_info "[dry-run] Would create directory: ${DEPTH_OUTPUT_DIR}"
        log_info "[dry-run] Would run: Rscript ${DEPTH_R_SCRIPT} ${SITE_METRICS_PATH} ${DEPTH_OUTPUT_DIR} ${PLOT_IMAGE_FORMAT}"
    else
        mkdir -p "${DEPTH_OUTPUT_DIR}"
        missing_depth_plots=()
        for chr in "${CHR_LABELS[@]}"; do
            expected_plot="${DEPTH_OUTPUT_DIR}/${chr}_depth_vs_position.${PLOT_IMAGE_FORMAT}"
            if [ ! -f "${expected_plot}" ]; then
                missing_depth_plots+=("${chr}")
            fi
        done
        if [ ${#missing_depth_plots[@]} -eq 0 ]; then
            log_info "All depth plots already exist; skipping generation."
        else
            log_info "Running depth plotting script..."
            log_info "Metrics: ${SITE_METRICS_PATH}"
            log_info "Output directory: ${DEPTH_OUTPUT_DIR}"
            if Rscript "${DEPTH_R_SCRIPT}" "${SITE_METRICS_PATH}" "${DEPTH_OUTPUT_DIR}" "${PLOT_IMAGE_FORMAT}"; then
                DEPTH_COUNT=$(find "${DEPTH_OUTPUT_DIR}" -maxdepth 1 -type f -name "*_depth_vs_position.${PLOT_IMAGE_FORMAT}" | awk 'END{print NR}')
                log_success "Generated ${DEPTH_COUNT} depth plots"
                log_info "Location: ${DEPTH_OUTPUT_DIR}/"
            else
                log_error "Depth plotting script failed"
                exit 1
            fi
        fi
    fi
else
    log_info "Skipping depth plots (mean depth metrics unavailable in this mode)."
fi

# =============================================================================
# STEP 4: GENERATE MISSINGNESS PLOTS
# =============================================================================

log_section "Step 4: Generate Missingness Plots"

MISSINGNESS_R_SCRIPT="${R_SCRIPTS_DIR}/plot_missingness_vs_position.R"

if ! check_file "${MISSINGNESS_R_SCRIPT}"; then
    log_error "R script not found: ${MISSINGNESS_R_SCRIPT}"
    exit 1
fi

if [ "${DRY_RUN}" = "true" ]; then
    log_info "[dry-run] Would create directory: ${MISSINGNESS_OUTPUT_DIR}"
    log_info "[dry-run] Would run: Rscript ${MISSINGNESS_R_SCRIPT} ${SITE_METRICS_PATH} ${MISSINGNESS_OUTPUT_DIR} ${PLOT_IMAGE_FORMAT}"
else
    mkdir -p "${MISSINGNESS_OUTPUT_DIR}"
    missing_missingness=()
    for chr in "${CHR_LABELS[@]}"; do
        expected_plot="${MISSINGNESS_OUTPUT_DIR}/${chr}_missingness_vs_position.${PLOT_IMAGE_FORMAT}"
        [ -f "${expected_plot}" ] || missing_missingness+=("${chr}")
    done
    if [ ${#missing_missingness[@]} -eq 0 ]; then
        log_info "All missingness plots already exist; skipping generation."
    else
        if Rscript "${MISSINGNESS_R_SCRIPT}" "${SITE_METRICS_PATH}" "${MISSINGNESS_OUTPUT_DIR}" "${PLOT_IMAGE_FORMAT}"; then
            MISSINGNESS_COUNT=$(find "${MISSINGNESS_OUTPUT_DIR}" -maxdepth 1 -type f -name "*_missingness_vs_position.${PLOT_IMAGE_FORMAT}" | awk 'END{print NR}')
            log_success "Generated ${MISSINGNESS_COUNT} missingness plots"
            log_info "Location: ${MISSINGNESS_OUTPUT_DIR}/"
        else
            log_error "Missingness plotting script failed"
            exit 1
        fi
    fi
fi

# =============================================================================
# STEP 5: GENERATE DEPTH VS MISSINGNESS SCATTER PLOTS
# =============================================================================

log_section "Step 5: Generate Depth vs Missingness Scatter Plots"

if [ "${HAS_DP}" = "true" ]; then
    DEPTH_MISS_R_SCRIPT="${R_SCRIPTS_DIR}/plot_depth_vs_missingness.R"

    if ! check_file "${DEPTH_MISS_R_SCRIPT}"; then
        log_error "R script not found: ${DEPTH_MISS_R_SCRIPT}"
        exit 1
    fi

    if [ "${DRY_RUN}" = "true" ]; then
        log_info "[dry-run] Would create directory: ${DEPTH_MISS_OUTPUT_DIR}"
        log_info "[dry-run] Would run: Rscript ${DEPTH_MISS_R_SCRIPT} ${SITE_METRICS_PATH} ${DEPTH_MISS_OUTPUT_DIR} ${PLOT_IMAGE_FORMAT}"
    else
        mkdir -p "${DEPTH_MISS_OUTPUT_DIR}"
        missing_depth_miss=()
        for chr in "${CHR_LABELS[@]}"; do
            expected_plot="${DEPTH_MISS_OUTPUT_DIR}/${chr}_depth_vs_missingness.${PLOT_IMAGE_FORMAT}"
            [ -f "${expected_plot}" ] || missing_depth_miss+=("${chr}")
        done
        if [ ${#missing_depth_miss[@]} -eq 0 ]; then
            log_info "All depth vs missingness plots already exist; skipping generation."
        else
            if Rscript "${DEPTH_MISS_R_SCRIPT}" "${SITE_METRICS_PATH}" "${DEPTH_MISS_OUTPUT_DIR}" "${PLOT_IMAGE_FORMAT}"; then
                DEPTH_MISS_COUNT=$(find "${DEPTH_MISS_OUTPUT_DIR}" -maxdepth 1 -type f -name "*_depth_vs_missingness.${PLOT_IMAGE_FORMAT}" | awk 'END{print NR}')
                log_success "Generated ${DEPTH_MISS_COUNT} depth vs missingness plots"
                log_info "Location: ${DEPTH_MISS_OUTPUT_DIR}/"
            else
                log_error "Depth vs missingness plotting script failed"
                exit 1
            fi
        fi
    fi
else
    log_info "Skipping depth vs missingness plots (mean depth metrics unavailable in this mode)."
fi

# =============================================================================
# STEP 6: GENERATE SITE QUALITY PLOTS
# =============================================================================

log_section "Step 6: Generate Site Quality Plots"

SITE_QUALITY_R_SCRIPT="${R_SCRIPTS_DIR}/plot_site_quality.R"

if ! check_file "${SITE_QUALITY_R_SCRIPT}"; then
    log_error "R script not found: ${SITE_QUALITY_R_SCRIPT}"
    exit 1
fi

if [ "${DRY_RUN}" = "true" ]; then
    log_info "[dry-run] Would create directory: ${SITE_QUALITY_OUTPUT_DIR}"
    log_info "[dry-run] Would run: Rscript ${SITE_QUALITY_R_SCRIPT} ${SITE_METRICS_PATH} ${SITE_QUALITY_OUTPUT_DIR} ${PLOT_IMAGE_FORMAT}"
else
    mkdir -p "${SITE_QUALITY_OUTPUT_DIR}"
    missing_quality=()
    for chr in "${CHR_LABELS[@]}"; do
        expected_plot="${SITE_QUALITY_OUTPUT_DIR}/${chr}_site_quality.${PLOT_IMAGE_FORMAT}"
        [ -f "${expected_plot}" ] || missing_quality+=("${chr}")
    done
    if [ ${#missing_quality[@]} -eq 0 ]; then
        log_info "All site quality plots already exist; skipping generation."
    else
        if Rscript "${SITE_QUALITY_R_SCRIPT}" "${SITE_METRICS_PATH}" "${SITE_QUALITY_OUTPUT_DIR}" "${PLOT_IMAGE_FORMAT}"; then
            SITE_QUALITY_COUNT=$(find "${SITE_QUALITY_OUTPUT_DIR}" -maxdepth 1 -type f -name "*_site_quality.${PLOT_IMAGE_FORMAT}" | awk 'END{print NR}')
            log_success "Generated ${SITE_QUALITY_COUNT} site quality plots"
            log_info "Location: ${SITE_QUALITY_OUTPUT_DIR}/"
        else
            log_error "Site quality plotting script failed"
            exit 1
        fi
    fi
fi

# =============================================================================
# STEP 7: GENERATE HETEROZYGOSITY PLOTS
# =============================================================================

log_section "Step 7: Generate Heterozygosity Plots"

HETEROZYGOSITY_R_SCRIPT="${R_SCRIPTS_DIR}/plot_heterozygosity.R"

if ! check_file "${HETEROZYGOSITY_R_SCRIPT}"; then
    log_error "R script not found: ${HETEROZYGOSITY_R_SCRIPT}"
    exit 1
fi

if [ "${DRY_RUN}" = "true" ]; then
    log_info "[dry-run] Would create directory: ${HETEROZYGOSITY_OUTPUT_DIR}"
    log_info "[dry-run] Would run: Rscript ${HETEROZYGOSITY_R_SCRIPT} ${SITE_METRICS_PATH} ${HETEROZYGOSITY_OUTPUT_DIR} ${PLOT_IMAGE_FORMAT}"
else
    mkdir -p "${HETEROZYGOSITY_OUTPUT_DIR}"
    missing_hetero=()
    for chr in "${CHR_LABELS[@]}"; do
        expected_plot="${HETEROZYGOSITY_OUTPUT_DIR}/${chr}_heterozygosity.${PLOT_IMAGE_FORMAT}"
        [ -f "${expected_plot}" ] || missing_hetero+=("${chr}")
    done
    if [ ${#missing_hetero[@]} -eq 0 ]; then
        log_info "All heterozygosity plots already exist; skipping generation."
    else
        if Rscript "${HETEROZYGOSITY_R_SCRIPT}" "${SITE_METRICS_PATH}" "${HETEROZYGOSITY_OUTPUT_DIR}" "${PLOT_IMAGE_FORMAT}"; then
            HETERO_COUNT=$(find "${HETEROZYGOSITY_OUTPUT_DIR}" -maxdepth 1 -type f -name "*_heterozygosity.${PLOT_IMAGE_FORMAT}" | awk 'END{print NR}')
            log_success "Generated ${HETERO_COUNT} heterozygosity plots"
            log_info "Location: ${HETEROZYGOSITY_OUTPUT_DIR}/"
        else
            log_error "Heterozygosity plotting script failed"
            exit 1
        fi
    fi
fi

# =============================================================================
# STEP 8: GENERATE QUALITY-BY-DEPTH PLOTS
# =============================================================================

log_section "Step 8: Generate Quality-by-Depth Plots"

if [ "${HAS_QD}" = "true" ]; then
    QD_R_SCRIPT="${R_SCRIPTS_DIR}/plot_quality_by_depth.R"

    if ! check_file "${QD_R_SCRIPT}"; then
        log_error "R script not found: ${QD_R_SCRIPT}"
        exit 1
    fi

    if [ "${DRY_RUN}" = "true" ]; then
        log_info "[dry-run] Would create directory: ${QD_OUTPUT_DIR}"
        log_info "[dry-run] Would run: Rscript ${QD_R_SCRIPT} ${SITE_METRICS_PATH} ${QD_OUTPUT_DIR} ${PLOT_IMAGE_FORMAT}"
    else
        mkdir -p "${QD_OUTPUT_DIR}"
        if [ -f "${QD_OUTPUT_FILE}" ]; then
            log_info "Quality-by-depth plot already exists; skipping generation."
        else
            if Rscript "${QD_R_SCRIPT}" "${SITE_METRICS_PATH}" "${QD_OUTPUT_DIR}" "${PLOT_IMAGE_FORMAT}"; then
                if [ -f "${QD_OUTPUT_FILE}" ]; then
                    log_success "Quality-by-depth plot created: ${QD_OUTPUT_FILE}"
                else
                    log_warning "Quality-by-depth plot script completed but output not found at ${QD_OUTPUT_FILE}"
                fi
            else
                log_error "Quality-by-depth plotting script failed"
                exit 1
            fi
        fi
    fi
else
    log_info "Skipping quality-by-depth histograms (QD metric unavailable in this mode)."
fi

# =============================================================================
# STEP 9: GENERATE CALL RATE HEAT MAPS
# =============================================================================

log_section "Step 9: Generate Call Rate Heat Maps"

CALL_RATE_R_SCRIPT="${R_SCRIPTS_DIR}/plot_call_rate_heatmap.R"

if ! check_file "${CALL_RATE_R_SCRIPT}"; then
    log_error "R script not found: ${CALL_RATE_R_SCRIPT}"
    exit 1
fi

if [ "${DRY_RUN}" = "true" ]; then
    log_info "[dry-run] Would create directory: ${CALL_RATE_OUTPUT_DIR}"
    log_info "[dry-run] Would run: Rscript ${CALL_RATE_R_SCRIPT} ${SITE_METRICS_PATH} ${CALL_RATE_OUTPUT_DIR} ${PLOT_IMAGE_FORMAT} ${CALL_RATE_HEATMAP_BINS}"
else
    mkdir -p "${CALL_RATE_OUTPUT_DIR}"
    missing_callrate=()
    for chr in "${CHR_LABELS[@]}"; do
        expected_plot="${CALL_RATE_OUTPUT_DIR}/${chr}_call_rate_heatmap.${PLOT_IMAGE_FORMAT}"
        [ -f "${expected_plot}" ] || missing_callrate+=("${chr}")
    done
    if [ ${#missing_callrate[@]} -eq 0 ]; then
        log_info "All call rate heat maps already exist; skipping generation."
    else
        if Rscript "${CALL_RATE_R_SCRIPT}" "${SITE_METRICS_PATH}" "${CALL_RATE_OUTPUT_DIR}" "${PLOT_IMAGE_FORMAT}" "${CALL_RATE_HEATMAP_BINS}"; then
            CALL_RATE_COUNT=$(find "${CALL_RATE_OUTPUT_DIR}" -maxdepth 1 -type f -name "*_call_rate_heatmap.${PLOT_IMAGE_FORMAT}" | awk 'END{print NR}')
            log_success "Generated ${CALL_RATE_COUNT} call rate heat maps"
            log_info "Location: ${CALL_RATE_OUTPUT_DIR}/"
        else
            log_error "Call rate heat map script failed"
            exit 1
        fi
    fi
fi

else
log_section "QC Stages Skipped"
log_info "QC mode not selected; skipping site metrics, depth/missingness/quality plots, AF distributions, and call rate heatmaps."
fi

# =============================================================================
# STEP 10: PERFORM PCA (PLINK2)
# =============================================================================

if [ "${RUN_DUP_CHECK}" = "true" ]; then
    log_section "Step 10: Duplicate Check (KING)"
    PCA_SCRIPT="${MODULE_DIR}/templates/plink2_PCA.sh"
    if ! check_file "${PCA_SCRIPT}"; then
        log_error "Duplicate check helper not found: ${PCA_SCRIPT}"
        exit 1
    fi

    mkdir -p "${PCA_OUTPUT_DIR}"
    if (cd "${PCA_OUTPUT_DIR}" && bash "${PCA_SCRIPT}" "${WORK_DIR}" "${R_SCRIPTS_DIR}" "${PLINK2_BIN_PATH}" "${BCFTOOLS_BIN_PATH}" "false" "false" "${LABEL_SIZE}" "false" "${PCA_MERGED_PATTERN}" "${DUPLICATE_MODE}" "${DUPLICATE_KING_THRESHOLD}" "duplicate"); then
        log_success "Duplicate check completed"
        log_info "Outputs located at ${PCA_OUTPUT_DIR}/ (king_duplicate_pairs.tsv, king_duplicate_samples.tsv)"
    else
        log_error "Duplicate check failed"
        exit 1
    fi
fi

# STEP 11: PERFORM PCA (PLINK2)
# =============================================================================

log_section "Step 11: Principal Component Analysis"

if [ "${RUN_PCA}" = "true" ]; then
    PCA_SCRIPT="${MODULE_DIR}/templates/plink2_PCA.sh"
    if ! check_file "${PCA_SCRIPT}"; then
        PCA_SKIP_REASON="Required PCA helper script missing"
        log_error "PCA script not found: ${PCA_SCRIPT}"
        exit 1
    fi

    if [ "${DRY_RUN}" = "true" ]; then
        PCA_SKIP_REASON="Dry-run mode (no PCA executed)"
        log_info "[dry-run] Would run PCA analysis in ${PCA_OUTPUT_DIR}"
    else
        if ! command -v "${PLINK2_BIN_PATH}" >/dev/null 2>&1; then
            PCA_SKIP_REASON="plink2 binary '${PLINK2_BIN_PATH}' not found in PATH"
            log_error "${PCA_SKIP_REASON}"
            exit 1
        fi
        if ! command -v "${BCFTOOLS_BIN_PATH}" >/dev/null 2>&1; then
            PCA_SKIP_REASON="bcftools binary '${BCFTOOLS_BIN_PATH}' not found in PATH"
            log_error "${PCA_SKIP_REASON}"
            exit 1
        fi

        mkdir -p "${PCA_OUTPUT_DIR}"
        if (cd "${PCA_OUTPUT_DIR}" && bash "${PCA_SCRIPT}" "${WORK_DIR}" "${R_SCRIPTS_DIR}" "${PLINK2_BIN_PATH}" "${BCFTOOLS_BIN_PATH}" "${REMOVE_RELATIVES}" "${SHOW_LABELS}" "${LABEL_SIZE}" "${USE_GGREPEL}" "${PCA_MERGED_PATTERN}" "${DUPLICATE_MODE}" "${DUPLICATE_KING_THRESHOLD}" "pca"); then
            PCA_EXECUTED="true"
            PCA_SKIP_REASON=""
            log_success "PCA analysis completed"
            log_info "PCA outputs located at ${PCA_OUTPUT_DIR}/"
        else
            PCA_SKIP_REASON="plink2_PCA.sh exited with a non-zero status"
            log_error "PCA analysis failed"
            exit 1
        fi
    fi
else
    if [ "${RUN_DUP_CHECK}" = "true" ]; then
        PCA_SKIP_REASON="Duplicate-check mode selected"
    elif [ "${RUN_QC}" = "true" ]; then
        PCA_SKIP_REASON="QC-only mode selected"
    else
        PCA_SKIP_REASON="PCA mode not selected"
    fi
    log_info "PCA analysis skipped (${PCA_SKIP_REASON})."
fi

# =============================================================================
# SUMMARY / DRY-RUN EXIT
# =============================================================================

if [ "${DRY_RUN}" = "true" ]; then
    log_section "Dry-Run Summary"
    echo "✅ Dry-run complete. No files were created."
    echo "   Planned outputs include:"
    echo "   • Site metrics TSV: ${SITE_METRICS_PATH}"
    if [ "${HAS_DP}" = "true" ]; then
        echo "   • Mean depth TSV:   ${MEAN_DP_PATH}"
    else
        echo "   • Mean depth TSV:   (skipped when depth metrics are unavailable)"
    fi
    echo "   • Plot directories: ${MISSINGNESS_OUTPUT_DIR}, ${SITE_QUALITY_OUTPUT_DIR}, ${HETEROZYGOSITY_OUTPUT_DIR}, ${CALL_RATE_OUTPUT_DIR}"
    if [ "${HAS_DP}" = "true" ]; then
        echo "                        ${DEPTH_OUTPUT_DIR}, ${DEPTH_MISS_OUTPUT_DIR}, ${QD_OUTPUT_DIR}"
    else
        echo "                        (depth, depth-vs-missingness, and QD plots skipped)"
    fi
    if [ "${RUN_PCA}" = "true" ]; then
        echo "   • PCA outputs:      ${PCA_OUTPUT_DIR}"
    fi
    echo ""
    echo "Review the log above to ensure paths and dependencies are correct."
    exit 0
fi

# =============================================================================
# SUMMARY
# =============================================================================

log_section "Pipeline Complete - Summary"

echo "📊 Analysis Results:"
echo ""
echo "1. Site Metrics TSV:"
if [ -f "${SITE_METRICS_PATH}" ]; then
    echo "   ✅ ${SITE_METRICS_PATH}"
    echo "      Lines: $(($(wc -l < "${SITE_METRICS_PATH}") - 1)) (excluding header)"
else
    if [ "${RUN_QC}" != "true" ]; then
        echo "   ⚠️  Not generated (QC mode not selected)."
    else
        echo "   ⚠️  Expected file not found at ${SITE_METRICS_PATH}"
    fi
fi
echo ""

echo "2. Filtered SNP VCFs:"
SNPS_COUNT=$(find "${WORK_DIR}" -maxdepth 1 -type f -name '*_snps.vcf.gz' | awk 'END{print NR}')
if [ "${SNPS_COUNT}" -eq "${EXPECTED_CHROM_COUNT}" ] && [ "${EXPECTED_CHROM_COUNT}" -gt 0 ]; then
    echo "   ✅ ${WORK_DIR}/"
    echo "      Files: ${SNPS_COUNT}/${EXPECTED_CHROM_COUNT} (_snps.vcf.gz)"
elif [ "${SNPS_COUNT}" -gt 0 ]; then
    echo "   ⚠️  ${SNPS_COUNT}/${EXPECTED_CHROM_COUNT} filtered SNP VCFs present"
else
    echo "   ⚠️  No *_snps.vcf.gz files detected in ${WORK_DIR}/"
fi
echo ""

echo "3. Per-Chromosome Metrics Directory:"
if [ -d "${METRICS_BY_CHROM_PATH}" ]; then
    FILE_COUNT=$(find "${METRICS_BY_CHROM_PATH}" -maxdepth 1 -type f -name '*_metrics.tsv' | awk 'END{print NR}')
    status_icon="✅"
    if [ "${FILE_COUNT}" -ne "${EXPECTED_CHROM_COUNT}" ]; then
        status_icon="⚠️"
    fi
    echo "   ${status_icon} ${METRICS_BY_CHROM_PATH}/"
    echo "      Files: ${FILE_COUNT}/${EXPECTED_CHROM_COUNT}"
else
    if [ "${RUN_QC}" != "true" ]; then
        echo "   ⚠️  Not generated (QC mode not selected)."
    else
        echo "   ⚠️  Metrics directory not found at ${METRICS_BY_CHROM_PATH}/"
    fi
fi
echo ""

echo "4. Mean Depth TSV:"
if [ "${HAS_DP}" = "true" ]; then
    if [ -f "${MEAN_DP_PATH}" ]; then
        echo "   ✅ ${MEAN_DP_PATH}"
        echo "      Lines: $(wc -l < "${MEAN_DP_PATH}")"
    else
        echo "   ⚠️  Expected file not found at ${MEAN_DP_PATH}"
    fi
else
    echo "   ⚠️  Skipped (mean depth metrics unavailable in this mode)"
fi
echo ""

echo "5. Mean Depth Plots:"
if [ "${HAS_DP}" = "true" ]; then
    if [ -d "${DEPTH_OUTPUT_DIR}" ]; then
        DEPTH_COUNT=$(find "${DEPTH_OUTPUT_DIR}" -maxdepth 1 -type f -name "*_depth_vs_position.${PLOT_IMAGE_FORMAT}" | awk 'END{print NR}')
        status_icon="✅"
        if [ "${DEPTH_COUNT}" -ne "${EXPECTED_CHROM_COUNT}" ]; then
            status_icon="⚠️"
        fi
        echo "   ${status_icon} ${DEPTH_OUTPUT_DIR}/"
        echo "      Plots: ${DEPTH_COUNT}/${EXPECTED_CHROM_COUNT}"
    else
        echo "   ⚠️  Depth plot directory not found at ${DEPTH_OUTPUT_DIR}/"
    fi
else
    echo "   ⚠️  Skipped (mean depth metrics unavailable in this mode)"
fi
echo ""

echo "6. Missingness Plots:"
if [ -d "${MISSINGNESS_OUTPUT_DIR}" ]; then
    MISS_COUNT=$(find "${MISSINGNESS_OUTPUT_DIR}" -maxdepth 1 -type f -name "*_missingness_vs_position.${PLOT_IMAGE_FORMAT}" | awk 'END{print NR}')
    status_icon="✅"
    if [ "${MISS_COUNT}" -ne "${EXPECTED_CHROM_COUNT}" ]; then
        status_icon="⚠️"
    fi
    echo "   ${status_icon} ${MISSINGNESS_OUTPUT_DIR}/"
    echo "      Plots: ${MISS_COUNT}/${EXPECTED_CHROM_COUNT}"
else
    if [ "${RUN_QC}" != "true" ]; then
        echo "   ⚠️  Not generated (QC mode not selected)."
    else
        echo "   ⚠️  Missingness plot directory not found at ${MISSINGNESS_OUTPUT_DIR}/"
    fi
fi
echo ""

echo "7. Allele Frequency Plots:"
if [ "${RUN_QC}" = "true" ]; then
    if [ -d "${AF_OUTPUT_DIR}" ]; then
        AF_COUNT=$(find "${AF_OUTPUT_DIR}" -maxdepth 1 -type f -name "*_af_distribution.${PLOT_IMAGE_FORMAT}" | awk 'END{print NR}')
        status_icon="✅"
        if [ "${AF_COUNT}" -lt "${EXPECTED_CHROM_COUNT}" ]; then
            status_icon="⚠️"
        fi
        echo "   ${status_icon} ${AF_OUTPUT_DIR}/"
        echo "      Plots: ${AF_COUNT} (including combined)"
    else
        echo "   ⚠️  Allele frequency directory not found at ${AF_OUTPUT_DIR}/"
    fi
else
    echo "   ⚠️  Skipped (QC mode not selected)"
fi
echo ""

echo "8. Depth vs Missingness Plots:"
if [ "${HAS_DP}" = "true" ]; then
    if [ -d "${DEPTH_MISS_OUTPUT_DIR}" ]; then
        DM_COUNT=$(find "${DEPTH_MISS_OUTPUT_DIR}" -maxdepth 1 -type f -name "*_depth_vs_missingness.${PLOT_IMAGE_FORMAT}" | awk 'END{print NR}')
        status_icon="✅"
        if [ "${DM_COUNT}" -ne "${EXPECTED_CHROM_COUNT}" ]; then
            status_icon="⚠️"
        fi
        echo "   ${status_icon} ${DEPTH_MISS_OUTPUT_DIR}/"
        echo "      Plots: ${DM_COUNT}/${EXPECTED_CHROM_COUNT}"
    else
        echo "   ⚠️  Depth-vs-missingness directory not found at ${DEPTH_MISS_OUTPUT_DIR}/"
    fi
else
    echo "   ⚠️  Skipped (mean depth metrics unavailable in this mode)"
fi
echo ""

echo "9. Site Quality Plots:"
if [ -d "${SITE_QUALITY_OUTPUT_DIR}" ]; then
    QUAL_COUNT=$(find "${SITE_QUALITY_OUTPUT_DIR}" -maxdepth 1 -type f -name "*_site_quality.${PLOT_IMAGE_FORMAT}" | awk 'END{print NR}')
    status_icon="✅"
    if [ "${QUAL_COUNT}" -ne "${EXPECTED_CHROM_COUNT}" ]; then
        status_icon="⚠️"
    fi
    echo "   ${status_icon} ${SITE_QUALITY_OUTPUT_DIR}/"
    echo "      Plots: ${QUAL_COUNT}/${EXPECTED_CHROM_COUNT}"
else
    if [ "${RUN_QC}" != "true" ]; then
        echo "   ⚠️  Not generated (QC mode not selected)."
    else
        echo "   ⚠️  Site quality directory not found at ${SITE_QUALITY_OUTPUT_DIR}/"
    fi
fi
echo ""

echo "10. Heterozygosity Plots:"
if [ -d "${HETEROZYGOSITY_OUTPUT_DIR}" ]; then
    HET_COUNT=$(find "${HETEROZYGOSITY_OUTPUT_DIR}" -maxdepth 1 -type f -name "*_heterozygosity.${PLOT_IMAGE_FORMAT}" | awk 'END{print NR}')
    status_icon="✅"
    if [ "${HET_COUNT}" -ne "${EXPECTED_CHROM_COUNT}" ]; then
        status_icon="⚠️"
    fi
    echo "   ${status_icon} ${HETEROZYGOSITY_OUTPUT_DIR}/"
    echo "      Plots: ${HET_COUNT}/${EXPECTED_CHROM_COUNT}"
else
    if [ "${RUN_QC}" != "true" ]; then
        echo "   ⚠️  Not generated (QC mode not selected)."
    else
        echo "   ⚠️  Heterozygosity directory not found at ${HETEROZYGOSITY_OUTPUT_DIR}/"
    fi
fi
echo ""

echo "11. Quality-by-Depth Plot:"
if [ "${HAS_QD}" = "true" ]; then
    if [ -d "${QD_OUTPUT_DIR}" ]; then
        if [ -f "${QD_OUTPUT_DIR}/quality_by_depth_hist.${PLOT_IMAGE_FORMAT}" ]; then
            echo "   ✅ ${QD_OUTPUT_DIR}/quality_by_depth_hist.${PLOT_IMAGE_FORMAT}"
        else
            echo "   ⚠️  Expected plot not found in ${QD_OUTPUT_DIR}/"
        fi
    else
        echo "   ⚠️  Quality-by-depth directory not found at ${QD_OUTPUT_DIR}/"
    fi
else
    echo "   ⚠️  Skipped (QD metric unavailable in this mode)"
fi
echo ""

echo "12. Call Rate Heat Maps:"
if [ -d "${CALL_RATE_OUTPUT_DIR}" ]; then
    HEAT_COUNT=$(find "${CALL_RATE_OUTPUT_DIR}" -maxdepth 1 -type f -name "*_call_rate_heatmap.${PLOT_IMAGE_FORMAT}" | awk 'END{print NR}')
    status_icon="✅"
    if [ "${HEAT_COUNT}" -ne "${EXPECTED_CHROM_COUNT}" ]; then
        status_icon="⚠️"
    fi
    echo "   ${status_icon} ${CALL_RATE_OUTPUT_DIR}/"
    echo "      Plots: ${HEAT_COUNT}/${EXPECTED_CHROM_COUNT}"
else
    if [ "${RUN_QC}" != "true" ]; then
        echo "   ⚠️  Not generated (QC mode not selected)."
    else
        echo "   ⚠️  Call-rate heat map directory not found at ${CALL_RATE_OUTPUT_DIR}/"
    fi
fi
echo ""

echo "13. PCA Outputs:"
if [ "${RUN_PCA}" = "true" ]; then
    if [ "${PCA_EXECUTED}" = "true" ]; then
        if [ -d "${PCA_OUTPUT_DIR}" ]; then
            PCA_FILE_COUNT=$(ls -1 "${PCA_OUTPUT_DIR}" 2>/dev/null | wc -l)
            echo "   ✅ ${PCA_OUTPUT_DIR}/"
            echo "      Files: ${PCA_FILE_COUNT}"
        else
            echo "   ⚠️  PCA reported success but directory is missing: ${PCA_OUTPUT_DIR}/"
        fi
    else
        echo "   ⚠️  PCA not executed: ${PCA_SKIP_REASON:-see logs for details}"
    fi
else
    echo "   ⚠️  Skipped (PCA analysis disabled)"
fi
echo ""

if [ "${RUN_DUP_CHECK}" = "true" ]; then
    echo "14. Duplicate Check Outputs:"
    if [ -d "${PCA_OUTPUT_DIR}" ]; then
        pairs="${PCA_OUTPUT_DIR}/king_duplicate_pairs.tsv"
        samples="${PCA_OUTPUT_DIR}/king_duplicate_samples.tsv"
        if [ -f "${pairs}" ]; then
            echo "   ✅ ${pairs}"
        else
            echo "   ⚠️  Duplicate pairs file not found at ${pairs}"
        fi
        if [ -f "${samples}" ]; then
            echo "   ✅ ${samples}"
        else
            echo "   ⚠️  Duplicate sample list not found at ${samples}"
        fi
    else
        echo "   ⚠️  Duplicate output directory not found at ${PCA_OUTPUT_DIR}/"
    fi
    echo ""
fi

log_success "All tasks completed successfully!"
log_info "Pipeline finished at $(date)"

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
