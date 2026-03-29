#!/bin/bash
################################################################################
# Master VCF Analysis Pipeline (v2 -- report-data)
# =============================================================================
# Compute-everything, decide-later workflow.
# The cluster extracts every metric, runs KING (all pairwise), and PCA on all
# samples, then exports a self-contained Parquet report package for local
# Quarto rendering.  No plots, no threshold filtering, no sample removal.
#
# Outputs:
#   STEP1D_CACHE_DIR  -- intermediate artefacts (SNP VCFs, metrics TSVs,
#                        combined VCF, PLINK files, KING .kin0)
#   STEP1D_PACKAGE_DIR -- final Parquet report package (tar'd)
#
# Usage:
#   Interactive: bash master_vcf_analysis.sh [--beagle] [--dry-run]
#   SLURM:       sbatch master_vcf_analysis.sh [--beagle] [--dry-run]
#
# NOTE: This script always filters inputs to biallelic SNPs
#       (bcftools view -m2 -M2 -v snps).
################################################################################

set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODULE_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
if [ -n "${PIPELINE_ROOT:-}" ]; then
    if [ -d "${PIPELINE_ROOT}" ]; then
        PIPELINE_ROOT="$(cd "${PIPELINE_ROOT}" && pwd)"
    else
        echo "[master_vcf_analysis] Warning: Provided PIPELINE_ROOT (${PIPELINE_ROOT}) does not exist; falling back to template-relative path." >&2
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
log_success() { log_info "$1"; }
log_warning() { log_warn "$1"; }

trap 'log_error "Unexpected failure in master_vcf_analysis.sh"; exit 1' ERR

usage() {
    cat <<'EOF'
Usage: master_vcf_analysis.sh [--dry-run] [--beagle] [--help]

Deprecated flags (accepted with warning, ignored):
  --qc, --PCA, --pca, --duplicate-check, --remove-relatives

Options:
  --dry-run, -n   Print planned actions without creating files
  --beagle        Treat input VCFs as Beagle-imputed (INFO tags: AF, DR2, IMP)
  --help, -h      Show this message and exit
EOF
}

DRY_RUN=false
BEAGLE_MODE=false

for arg in "$@"; do
    case "$arg" in
        --dry-run|-n)
            DRY_RUN=true
            ;;
        --beagle)
            BEAGLE_MODE=true
            ;;
        --qc|--PCA|--pca|--duplicate-check|--remove-relatives)
            log_warn "Flag '${arg}' is deprecated and ignored. Step1D now runs a combined report-data workflow."
            log_warn "This flag will be removed in a future version."
            ;;
        --help|-h)
            usage
            exit 0
            ;;
        -*)
            echo "Unknown option: ${arg}" >&2
            usage
            exit 1
            ;;
        *)
            echo "Unexpected argument: ${arg}" >&2
            echo "" >&2
            echo "This script does not accept directory arguments directly." >&2
            echo "Use environment variables instead:" >&2
            echo "" >&2
            echo "  export VCF_DIR=${arg}" >&2
            echo "  bash master_vcf_analysis.sh" >&2
            echo "" >&2
            usage
            exit 1
            ;;
    esac
done

# =============================================================================
# CONFIGURATION
# =============================================================================

VCF_DIR="${VCF_DIR:-$PWD}"
VCF_PATTERN="${VCF_PATTERN:-Chr%02d.vcf.gz}"
WORK_DIR="${WORK_DIR:-${VCF_DIR}}"
R_SCRIPTS_DIR="${R_SCRIPTS_DIR:-${MODULE_DIR}/Rscripts}"

DATASET_NAME="${DATASET_NAME:-$(basename "${VCF_DIR}")}"
STEP1D_CACHE_DIR="${STEP1D_CACHE_DIR:-${WORK_DIR}/${DATASET_NAME}_step1d_cache}"
STEP1D_PACKAGE_DIR="${STEP1D_PACKAGE_DIR:-${WORK_DIR}/${DATASET_NAME}_step1d_report_package}"
STEP1D_PACKAGE_TAR="${STEP1D_PACKAGE_TAR:-${STEP1D_PACKAGE_DIR}.tar.gz}"
STEP1D_PARQUET_COMPRESSION="${STEP1D_PARQUET_COMPRESSION:-snappy}"

SITE_METRICS_TSV="${SITE_METRICS_TSV:-variant_site_metrics.tsv}"
METRICS_BY_CHROM_DIR="${METRICS_BY_CHROM_DIR:-site_metrics_per_chromosome}"
SITE_METRICS_PATH="${STEP1D_CACHE_DIR}/${SITE_METRICS_TSV}"
METRICS_BY_CHROM_PATH="${STEP1D_CACHE_DIR}/${METRICS_BY_CHROM_DIR}"
MEAN_DP_TSV="${MEAN_DP_TSV:-SNP_site_meanDP.tsv}"
MEAN_DP_PATH="${STEP1D_CACHE_DIR}/${MEAN_DP_TSV}"

PCA_DIR_NAME="${STEP1D_PCA_DIR:-pca_analysis}"
PCA_OUTPUT_DIR="${STEP1D_CACHE_DIR}/${PCA_DIR_NAME}"

PLINK2_BIN_PATH="${PLINK2_BIN:-plink2}"
BCFTOOLS_BIN_PATH="${BCFTOOLS_BIN:-bcftools}"
PCA_MERGED_PATTERN="${STEP1D_PCA_MERGED_PATTERN:-*merged*.vcf.gz,*merge*.vcf.gz}"
export STEP1D_PCA_FORCE_CONCAT="${STEP1D_PCA_FORCE_CONCAT:-false}"
export STEP1D_PCA_MERGED_EXCLUDE_CHR="${STEP1D_PCA_MERGED_EXCLUDE_CHR:-true}"

CONDA_ENV="${CONDA_ENV:-rplot}"

HAS_DP=false
HAS_QD=false

# Build VCF target list
declare -a VCF_TARGETS=()
if [ -n "${VCF_INCLUDE_FILENAMES:-}" ]; then
    read -r -a VCF_TARGETS <<< "${VCF_INCLUDE_FILENAMES}"
else
    for chr_num in $(seq 0 17); do
        VCF_TARGETS+=("$(printf "${VCF_PATTERN}" "${chr_num}")")
    done
fi

declare -a EXISTING_VCF_TARGETS=()
for vcf_basename in "${VCF_TARGETS[@]}"; do
    if [ -f "${VCF_DIR}/${vcf_basename}" ]; then
        EXISTING_VCF_TARGETS+=("${vcf_basename}")
    fi
done
VCF_TARGETS=("${EXISTING_VCF_TARGETS[@]}")

if [ ${#VCF_TARGETS[@]} -eq 0 ]; then
    declare -a AUTO_VCFS=()
    while IFS= read -r -d '' candidate; do
        AUTO_VCFS+=("$(basename "${candidate}")")
    done < <(find "${VCF_DIR}" -maxdepth 1 -type f \( -name "*.vcf.gz" -o -name "*.vcf" \) -print0 | LC_ALL=C sort -z)

    if [ ${#AUTO_VCFS[@]} -eq 0 ]; then
        log_error "No VCF files found in ${VCF_DIR} (looked for *.vcf.gz and *.vcf)."
        exit 1
    fi

    VCF_TARGETS=("${AUTO_VCFS[@]}")
    log_warn "No Chr-pattern VCFs found (VCF_PATTERN=${VCF_PATTERN}). Auto-detected ${#VCF_TARGETS[@]} VCF(s)."
    log_info "Auto-detected VCFs: ${VCF_TARGETS[*]}"
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

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

log_section() {
    echo ""
    echo "================================================================="
    echo "$1"
    echo "================================================================="
    echo ""
}

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
        for pkg in "${missing[@]}"; do
            log_error "    conda install -n ${CONDA_ENV} r-${pkg}"
        done
        log_error "Or activate ${CONDA_ENV} and run:"
        log_error "    Rscript -e \"install.packages(c('${missing[*]}'), repos='https://cloud.r-project.org')\""
        exit 1
    fi
}

# =============================================================================
# SETUP
# =============================================================================

log_section "Step1D Report-Data Pipeline - Starting"

log_info "Pipeline started"
log_info "VCF Directory (read-only input): ${VCF_DIR}"
log_info "Cache Directory: ${STEP1D_CACHE_DIR}"
log_info "Package Directory: ${STEP1D_PACKAGE_DIR}"
log_info "R Scripts Directory: ${R_SCRIPTS_DIR}"
if [ "${BEAGLE_MODE}" = "true" ]; then
    log_info "Beagle mode enabled: using AF/DR2/IMP metrics and skipping depth-dependent columns."
fi

if [ "${DRY_RUN}" = "true" ]; then
    log_info "Dry-run mode enabled. No files will be created."
fi

cd "${WORK_DIR}" || {
    log_error "Cannot change to working directory: ${WORK_DIR}"
    exit 1
}

mkdir -p "${STEP1D_CACHE_DIR}"

# Load HPC modules one at a time, capturing absolute binary paths to avoid
# GCC toolchain conflicts between plink (gcc-11) and bcftools (gcc-12).
log_info "Loading required modules..."

if command -v module >/dev/null 2>&1; then
    module purge 2>/dev/null

    module load miniforge/25.3.0-3 2>/dev/null && log_success "Loaded miniforge module" \
        || log_warn "miniforge module not loaded (may already be available)"

    PLINK_MODULE="${PLINK_MODULE:-plink/2.00a3.6-gcc-11.3.0}"
    if module load "${PLINK_MODULE}" 2>/dev/null; then
        PLINK2_BIN_PATH="$(command -v plink2)"
        log_success "Loaded ${PLINK_MODULE} -> ${PLINK2_BIN_PATH}"
        module unload "${PLINK_MODULE}" 2>/dev/null
    fi

    BCFTOOLS_MODULE="${BCFTOOLS_MODULE:-bcftools/1.18-gcc-12.3.0}"
    if module load "${BCFTOOLS_MODULE}" 2>/dev/null; then
        BCFTOOLS_BIN_PATH="$(command -v bcftools)"
        log_success "Loaded ${BCFTOOLS_MODULE} -> ${BCFTOOLS_BIN_PATH}"
        module unload "${BCFTOOLS_MODULE}" 2>/dev/null
    fi
else
    log_warn "module command not available; assuming tools are already on PATH"
fi

# Verify required binaries are available
if ! [ -x "${PLINK2_BIN_PATH}" ] && ! command -v "${PLINK2_BIN_PATH}" >/dev/null 2>&1; then
    log_error "plink2 binary not found (tried: ${PLINK2_BIN_PATH})"
    exit 1
fi
if ! [ -x "${BCFTOOLS_BIN_PATH}" ] && ! command -v "${BCFTOOLS_BIN_PATH}" >/dev/null 2>&1; then
    log_error "bcftools binary not found (tried: ${BCFTOOLS_BIN_PATH})"
    exit 1
fi
log_info "plink2: ${PLINK2_BIN_PATH}"
log_info "bcftools: ${BCFTOOLS_BIN_PATH}"

# Add binary directories to PATH so child scripts (e.g. prepare_combined_for_pca.sh) find them
export PATH="$(dirname "${PLINK2_BIN_PATH}"):$(dirname "${BCFTOOLS_BIN_PATH}"):${PATH}"

# Activate conda for R
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

log_info "Checking required commands..."
check_command Rscript || exit 1
ensure_r_packages data.table arrow jsonlite
log_success "Required commands and packages available"

# =============================================================================
# STEP 1: BUILD SITE-LEVEL METRICS DATASET
# =============================================================================

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
        log_warning "Existing site metrics TSV appears malformed. Regenerating."
        REGENERATE_METRICS=true
    else
        if [ "${BEAGLE_MODE}" = "true" ]; then
            for marker in AF DR2 IMP; do
                if ! printf '%s\n' "${header_line}" | tr '\t' '\n' | grep -Fxq "${marker}"; then
                    log_warning "Existing site metrics missing '${marker}' column. Regenerating."
                    REGENERATE_METRICS=true
                    break
                fi
            done
            if [ "${REGENERATE_METRICS}" = false ] && printf '%s\n' "${header_line}" | tr '\t' '\n' | grep -Fxq "MEAN_DEPTH"; then
                log_warning "Existing site metrics include depth columns but --beagle mode requested. Regenerating."
                REGENERATE_METRICS=true
            fi
        else
            for marker in MEAN_DEPTH QD MQ; do
                if ! printf '%s\n' "${header_line}" | tr '\t' '\n' | grep -Fxq "${marker}"; then
                    log_warning "Existing site metrics missing '${marker}' column. Regenerating."
                    REGENERATE_METRICS=true
                    break
                fi
            done
            if [ "${REGENERATE_METRICS}" = false ] && printf '%s\n' "${header_line}" | tr '\t' '\n' | grep -Fxq "DR2"; then
                log_warning "Existing site metrics appear Beagle-mode. Regenerating for standard dataset."
                REGENERATE_METRICS=true
            fi
        fi
    fi
fi

if [ "${REGENERATE_METRICS}" = true ]; then
    rm -f "${SITE_METRICS_PATH}"
    find "${METRICS_BY_CHROM_PATH}" -maxdepth 1 -type f -name '*_metrics.tsv' -delete 2>/dev/null || true
    find "${STEP1D_CACHE_DIR}" -maxdepth 1 -type f -name '*_snps.vcf.gz*' -delete 2>/dev/null || true
fi

if [ -f "${SITE_METRICS_PATH}" ]; then
    log_success "Site metrics TSV already exists: ${SITE_METRICS_PATH}"
    log_info "Skipping regeneration of site metrics"
else
    if [ "${DRY_RUN}" = "true" ]; then
        log_info "[dry-run] Target VCF files: ${VCF_TARGETS[*]}"
        log_info "[dry-run] Would create filtered SNP VCFs at ${STEP1D_CACHE_DIR}/<chromosome>_snps.vcf.gz"
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
            filtered_output="${STEP1D_CACHE_DIR}/${chr_name}_snps.vcf.gz"
            metrics_output="${METRICS_BY_CHROM_PATH}/${chr_name}_metrics.tsv"
            tmp_metrics="${TEMP_DIR}/${chr_name}_metrics.tsv"

            if [ ! -f "${vcf_file}" ]; then
                log_warning "${chr_name}: VCF not found (${vcf_file}), skipping"
                continue
            fi

            if [ -f "${metrics_output}" ] && [ "${REGENERATE_METRICS}" = false ]; then
                log_info "${chr_name}: metrics TSV already present; skipping."
                METRIC_FILES+=("${metrics_output}")
                continue
            fi

            log_info "${chr_name}: filtering to biallelic SNPs"
            if ! "${BCFTOOLS_BIN_PATH}" view -m2 -M2 -v snps -Oz -o "${filtered_temp}" "${vcf_file}"; then
                log_error "${chr_name}: failed to filter biallelic SNPs"
                rm -rf "${TEMP_DIR}"
                exit 1
            fi
            if ! "${BCFTOOLS_BIN_PATH}" index -t "${filtered_temp}"; then
                log_error "${chr_name}: failed to index filtered SNP VCF"
                rm -rf "${TEMP_DIR}"
                exit 1
            fi
            mv -f "${filtered_temp}" "${filtered_output}"
            mv -f "${filtered_temp}.tbi" "${filtered_output}.tbi"
            log_info "${chr_name}: filtered SNP VCF -> ${filtered_output}"

            log_info "${chr_name}: extracting site metrics..."

            if [ "${BEAGLE_MODE}" = "true" ]; then
                query_format='%CHROM\t%POS\t%QUAL\t%INFO/AF\t%INFO/DR2\t%INFO/IMP\t[%GT\t]\n'
                if ! "${BCFTOOLS_BIN_PATH}" query -f "${query_format}" "${filtered_output}" | \
                    awk -v OFS='\t' -v header="${SITE_METRICS_HEADER}" '
                        BEGIN { print header; }
                        {
                            sub(/\t$/, "", $0);
                            chrom = $1; pos = $2; qual = $3;
                            af = $4; dr2 = $5; imp = $6;
                            total_genotypes = NF - 6;
                            called = 0; missing = 0; hetero = 0;
                            for (i = 7; i <= NF; i++) {
                                gt = $i;
                                if (gt == "" || gt == "." || gt == "./." || gt == ".|.") { missing++; }
                                else {
                                    called++;
                                    split(gt, alleles, /[\/|]/);
                                    if (length(alleles) >= 2) {
                                        if (alleles[1] != "." && alleles[2] != "." && alleles[1] != alleles[2]) hetero++;
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
                if ! "${BCFTOOLS_BIN_PATH}" query \
                    -f '%CHROM\t%POS\t%QUAL\t%INFO/QD\t%INFO/AC\t%INFO/AF\t%INFO/InbreedingCoeff\t%INFO/ExcessHet\t%INFO/MQ\t[%GT:%DP\t]\n' \
                    "${filtered_output}" | \
                    awk -v OFS='\t' -v header="${SITE_METRICS_HEADER}" '
                        BEGIN { print header; }
                        {
                            sub(/\t$/, "", $0);
                            chrom = $1; pos = $2; qual = $3;
                            qd = $4; ac = $5; af = $6;
                            inbreeding = $7; excesshet = $8; mq = $9;
                            sum_dp = 0; dp_non_missing = 0;
                            total = 0; missing = 0; hetero = 0; called = 0;
                            qd_value = (qd == "." || qd == "" || qd == "NA") ? "NA" : qd;
                            ac_value = (ac == "" || ac == "NA" || ac == ".") ? "NA" : ac;
                            af_value = (af == "" || af == "NA" || af == ".") ? "NA" : af;
                            inbreed_value = (inbreeding == "" || inbreeding == "NA" || inbreeding == ".") ? "NA" : inbreeding;
                            excesshet_value = (excesshet == "" || excesshet == "NA" || excesshet == ".") ? "NA" : excesshet;
                            mq_value = (mq == "" || mq == "NA" || mq == ".") ? "NA" : mq;
                            for (i = 10; i <= NF; i++) {
                                field = $i;
                                if (field == "") continue;
                                total++;
                                split(field, parts, ":");
                                gt = parts[1];
                                dp = (length(parts) > 1 ? parts[2] : "");
                                if (dp != "" && dp != "." && dp != "./." && dp != ".|." && dp ~ /^-?[0-9]+(\.[0-9]+)?$/) {
                                    sum_dp += dp; dp_non_missing++;
                                }
                                if (gt != "" && gt != "." && gt != "./." && gt != ".|.") {
                                    split(gt, alleles, /[\/|]/);
                                    if (length(alleles) >= 2) {
                                        if (alleles[1] != "." && alleles[2] != ".") {
                                            called++;
                                            if (alleles[1] != alleles[2]) hetero++;
                                        }
                                    }
                                } else { missing++; }
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
            log_error "No per-chromosome metrics files are available. Aborting."
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

# Derive mean depth TSV (non-Beagle only)
if [ "${HAS_DP}" = "true" ] && [ "${DRY_RUN}" != "true" ]; then
    if [ -f "${MEAN_DP_PATH}" ]; then
        log_success "Mean depth TSV available: ${MEAN_DP_PATH}"
    elif [ -f "${SITE_METRICS_PATH}" ]; then
        log_info "Deriving mean depth TSV from site metrics..."
        {
            printf "CHROM\tPOS\tMEAN_DEPTH\n"
            awk 'NR == 1 { next } { printf "%s\t%s\t%s\n", $1, $2, $10 }' "${SITE_METRICS_PATH}"
        } > "${MEAN_DP_PATH}"
        log_success "Mean depth TSV created: ${MEAN_DP_PATH}"
    fi
fi

# =============================================================================
# STEP 2: PREPARE COMBINED VCF FOR PCA
# =============================================================================

log_section "Step 2: Prepare Combined VCF for PCA"

PREPARE_COMBINED_SCRIPT="${PIPELINE_ROOT}/modules/step1d/bin/prepare_combined_for_pca.sh"
COMBINED_FOR_PCA="${STEP1D_CACHE_DIR}/combined_for_pca.vcf.gz"
COMBINED_STATS="${STEP1D_CACHE_DIR}/combined_for_pca.stats.txt"

if [ "${DRY_RUN}" = "true" ]; then
    log_info "[dry-run] Would prepare combined VCF for PCA at ${COMBINED_FOR_PCA}"
else
    if ! check_file "${PREPARE_COMBINED_SCRIPT}"; then
        log_error "prepare_combined_for_pca.sh not found; cannot continue."
        exit 1
    fi

    # If combined VCF exists in VCF_DIR but not in cache, move it
    if [ ! -f "${COMBINED_FOR_PCA}" ] && [ -f "${VCF_DIR}/combined_for_pca.vcf.gz" ] && [ "${VCF_DIR}" != "${STEP1D_CACHE_DIR}" ]; then
        log_info "Found combined VCF in VCF_DIR; moving to cache..."
        mv -f "${VCF_DIR}/combined_for_pca.vcf.gz" "${COMBINED_FOR_PCA}"
        [ -f "${VCF_DIR}/combined_for_pca.vcf.gz.csi" ] && mv -f "${VCF_DIR}/combined_for_pca.vcf.gz.csi" "${COMBINED_FOR_PCA}.csi"
        [ -f "${VCF_DIR}/combined_for_pca.vcf.gz.tbi" ] && mv -f "${VCF_DIR}/combined_for_pca.vcf.gz.tbi" "${COMBINED_FOR_PCA}.tbi"
        [ -f "${VCF_DIR}/combined_for_pca.stats.txt" ] && mv -f "${VCF_DIR}/combined_for_pca.stats.txt" "${COMBINED_STATS}"
    fi

    NEED_COMBINED="false"
    if [ ! -f "${COMBINED_FOR_PCA}" ]; then
        log_info "Combined VCF for PCA not found; will prepare it now."
        NEED_COMBINED="true"
    else
        for vcf_basename in "${VCF_TARGETS[@]}"; do
            vcf_path="${VCF_DIR}/${vcf_basename}"
            if [ -f "${vcf_path}" ] && [ "${vcf_path}" -nt "${COMBINED_FOR_PCA}" ]; then
                log_info "Source VCF ${vcf_basename} is newer than combined VCF; will regenerate."
                NEED_COMBINED="true"
                break
            fi
        done
    fi

    if [ "${NEED_COMBINED}" = "true" ]; then
        log_info "Preparing combined VCF for PCA..."
        if bash "${PREPARE_COMBINED_SCRIPT}" "${VCF_DIR}"; then
            # Script writes to VCF_DIR; move to cache if needed
            if [ -f "${VCF_DIR}/combined_for_pca.vcf.gz" ] && [ "${VCF_DIR}/combined_for_pca.vcf.gz" != "${COMBINED_FOR_PCA}" ]; then
                mv -f "${VCF_DIR}/combined_for_pca.vcf.gz" "${COMBINED_FOR_PCA}"
                [ -f "${VCF_DIR}/combined_for_pca.vcf.gz.csi" ] && mv -f "${VCF_DIR}/combined_for_pca.vcf.gz.csi" "${COMBINED_FOR_PCA}.csi"
                [ -f "${VCF_DIR}/combined_for_pca.vcf.gz.tbi" ] && mv -f "${VCF_DIR}/combined_for_pca.vcf.gz.tbi" "${COMBINED_FOR_PCA}.tbi"
                [ -f "${VCF_DIR}/combined_for_pca.stats.txt" ] && mv -f "${VCF_DIR}/combined_for_pca.stats.txt" "${COMBINED_STATS}"
            fi
            log_success "Combined VCF prepared: ${COMBINED_FOR_PCA}"
        else
            log_error "Failed to prepare combined VCF for PCA"
            exit 1
        fi
    else
        log_info "Combined VCF already exists and is up-to-date: ${COMBINED_FOR_PCA}"
        if [ ! -f "${COMBINED_STATS}" ]; then
            log_info "Stats file missing; generating bcftools stats..."
            if "${BCFTOOLS_BIN_PATH}" stats "${COMBINED_FOR_PCA}" > "${COMBINED_STATS}"; then
                log_success "Stats generated for existing combined VCF"
            else
                log_warn "Failed to generate stats; continuing without stats summary"
            fi
        fi
    fi

    if [ -f "${COMBINED_STATS}" ]; then
        snp_count=$(grep "^SN" "${COMBINED_STATS}" | grep "number of SNPs:" | awk '{print $NF}' || echo "unknown")
        log_info "Total SNPs in combined VCF: ${snp_count}"
    fi
fi

# =============================================================================
# STEP 3: RUN PLINK2 PCA PIPELINE (import, QC, KING, LD prune, PCA)
# =============================================================================

log_section "Step 3: PLINK2 PCA Pipeline"

PCA_SCRIPT="${PIPELINE_ROOT}/modules/step1d/templates/plink2_PCA.sh"
if ! check_file "${PCA_SCRIPT}"; then
    log_error "PCA helper script not found: ${PCA_SCRIPT}"
    exit 1
fi

if [ "${DRY_RUN}" = "true" ]; then
    log_info "[dry-run] Would run PCA pipeline in ${PCA_OUTPUT_DIR}"
else
    mkdir -p "${PCA_OUTPUT_DIR}"
    if (cd "${PCA_OUTPUT_DIR}" && bash "${PCA_SCRIPT}" "${STEP1D_CACHE_DIR}" "${R_SCRIPTS_DIR}" "${PLINK2_BIN_PATH}" "${BCFTOOLS_BIN_PATH}" "${PCA_OUTPUT_DIR}"); then
        log_success "PCA pipeline completed"
        log_info "PCA outputs: ${PCA_OUTPUT_DIR}/"
    else
        log_error "PCA pipeline failed"
        exit 1
    fi
fi

# =============================================================================
# STEP 4: EXPORT PARQUET REPORT PACKAGE
# =============================================================================

log_section "Step 4: Export Parquet Report Package"

EXPORT_SCRIPT="${R_SCRIPTS_DIR}/export_parquet_package.R"

if [ "${DRY_RUN}" = "true" ]; then
    log_info "[dry-run] Would export Parquet package to ${STEP1D_PACKAGE_DIR}"
    log_info "[dry-run] Would create tarball: ${STEP1D_PACKAGE_TAR}"
else
    if ! check_file "${EXPORT_SCRIPT}"; then
        log_error "export_parquet_package.R not found at ${EXPORT_SCRIPT}"
        exit 1
    fi

    BEAGLE_FLAG_VALUE="false"
    if [ "${BEAGLE_MODE}" = "true" ]; then
        BEAGLE_FLAG_VALUE="true"
    fi

    log_info "Running Parquet exporter..."
    if Rscript "${EXPORT_SCRIPT}" \
        --cache-dir "${STEP1D_CACHE_DIR}" \
        --package-dir "${STEP1D_PACKAGE_DIR}" \
        --dataset-name "${DATASET_NAME}" \
        --beagle "${BEAGLE_FLAG_VALUE}" \
        --compression "${STEP1D_PARQUET_COMPRESSION}"; then
        log_success "Parquet package exported: ${STEP1D_PACKAGE_DIR}"
    else
        log_error "Parquet export failed"
        exit 1
    fi

    log_info "Creating tarball..."
    tar -czf "${STEP1D_PACKAGE_TAR}" -C "$(dirname "${STEP1D_PACKAGE_DIR}")" "$(basename "${STEP1D_PACKAGE_DIR}")"
    log_success "Report package tarball: ${STEP1D_PACKAGE_TAR}"

    TAR_SIZE=$(du -h "${STEP1D_PACKAGE_TAR}" | cut -f1)
    log_info "Tarball size: ${TAR_SIZE}"
fi

# =============================================================================
# SUMMARY
# =============================================================================

if [ "${DRY_RUN}" = "true" ]; then
    log_section "Dry-Run Summary"
    echo "Dry-run complete. No files were created."
    echo "   Planned outputs include:"
    echo "   - Cache dir:       ${STEP1D_CACHE_DIR}"
    echo "   - Site metrics:    ${SITE_METRICS_PATH}"
    echo "   - Combined VCF:    ${COMBINED_FOR_PCA}"
    echo "   - PCA outputs:     ${PCA_OUTPUT_DIR}"
    echo "   - Report package:  ${STEP1D_PACKAGE_DIR}"
    echo "   - Tarball:         ${STEP1D_PACKAGE_TAR}"
    exit 0
fi

log_section "Pipeline Complete - Summary"

echo "Report-Data Pipeline Results:"
echo ""
echo "1. Site Metrics TSV:"
if [ -f "${SITE_METRICS_PATH}" ]; then
    echo "   OK ${SITE_METRICS_PATH}"
    echo "      Lines: $(($(wc -l < "${SITE_METRICS_PATH}") - 1)) (excluding header)"
else
    echo "   WARN Expected file not found at ${SITE_METRICS_PATH}"
fi
echo ""

echo "2. Per-Chromosome Metrics:"
if [ -d "${METRICS_BY_CHROM_PATH}" ]; then
    FILE_COUNT=$(find "${METRICS_BY_CHROM_PATH}" -maxdepth 1 -type f -name '*_metrics.tsv' | awk 'END{print NR}')
    echo "   OK ${METRICS_BY_CHROM_PATH}/ (${FILE_COUNT}/${EXPECTED_CHROM_COUNT} files)"
else
    echo "   WARN Metrics directory not found"
fi
echo ""

echo "3. Combined VCF:"
if [ -f "${COMBINED_FOR_PCA}" ]; then
    echo "   OK ${COMBINED_FOR_PCA}"
else
    echo "   WARN ${COMBINED_FOR_PCA} not found"
fi
echo ""

echo "4. PCA Outputs:"
if [ -d "${PCA_OUTPUT_DIR}" ]; then
    PCA_FILE_COUNT=$(ls -1 "${PCA_OUTPUT_DIR}" 2>/dev/null | wc -l)
    echo "   OK ${PCA_OUTPUT_DIR}/ (${PCA_FILE_COUNT} files)"
else
    echo "   WARN PCA directory not found"
fi
echo ""

echo "5. Report Package:"
if [ -f "${STEP1D_PACKAGE_TAR}" ]; then
    TAR_SIZE=$(du -h "${STEP1D_PACKAGE_TAR}" | cut -f1)
    echo "   OK ${STEP1D_PACKAGE_TAR} (${TAR_SIZE})"
else
    echo "   WARN Tarball not found at ${STEP1D_PACKAGE_TAR}"
fi
echo ""

log_success "All tasks completed successfully!"
log_info "Pipeline finished at $(date)"

echo ""
echo "================================================================="
