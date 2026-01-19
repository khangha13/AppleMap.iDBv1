#!/bin/bash
################################################################################
# SNPs-per-sample count table from per-chromosome VCFs
#
# Generates a single gzip-compressed CSV with per-VCF (CHR) and per-sample TOTAL
# rows for called biallelic PASS SNP genotypes.
#
# Usage:
#   Interactive: bash snp_counts_per_sample_from_chr_vcfs.sh --vcf-dir /path/to/vcfs
#   SLURM: sbatch --export=ALL,VCF_DIR=/path/to/vcfs snp_counts_per_sample_from_chr_vcfs.sh
#
# R quick start:
#   dt <- data.table::fread("snps_per_sample_counts.csv.gz")
#   dt_total <- dt[chr_label=="TOTAL"]
#   ggplot2::ggplot(dt_total, ggplot2::aes(y=n_called, x="")) +
#     ggplot2::geom_boxplot() +
#     ggplot2::geom_jitter(width=0.15, alpha=0.3)
################################################################################

set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_ROOT_FALLBACK="$(cd "${SCRIPT_DIR}/.." && pwd)"
if [ -n "${PIPELINE_ROOT:-}" ]; then
    if [ -d "${PIPELINE_ROOT}" ]; then
        PIPELINE_ROOT="$(cd "${PIPELINE_ROOT}" && pwd)"
    else
        echo "[WARN] Provided PIPELINE_ROOT (${PIPELINE_ROOT}) does not exist; using script-relative path." >&2
        PIPELINE_ROOT="${PIPELINE_ROOT_FALLBACK}"
    fi
else
    PIPELINE_ROOT="${PIPELINE_ROOT_FALLBACK}"
fi
export PIPELINE_ROOT

if [ -f "${PIPELINE_ROOT}/config/pipeline_config.sh" ]; then
    # shellcheck disable=SC1091
    source "${PIPELINE_ROOT}/config/pipeline_config.sh"
fi
if [ -f "${PIPELINE_ROOT}/lib/logging.sh" ]; then
    # shellcheck disable=SC1091
    source "${PIPELINE_ROOT}/lib/logging.sh"
    init_logging "snp_counts_per_sample" "pipeline" "${DATASET_NAME:-${DATASET:-}}"
fi

if ! declare -F log_info >/dev/null 2>&1; then
    log_info() { echo "[INFO] $*" >&2; }
    log_warn() { echo "[WARN] $*" >&2; }
    log_error() { echo "[ERROR] $*" >&2; }
fi

trap 'log_error "Unexpected failure in snp_counts_per_sample_from_chr_vcfs.sh"; exit 1' ERR

usage() {
    cat <<'EOF_USAGE'
Usage: snp_counts_per_sample_from_chr_vcfs.sh [options]

Options:
  --vcf-dir <dir>        Directory to scan for *.vcf.gz (env: VCF_DIR, default: PWD)
  --work-dir <dir>       Output directory (env: WORK_DIR, default: VCF_DIR)
  --out <path>           Output path (env: OUT_FILE, default: WORK_DIR/snps_per_sample_counts.csv.gz)
  --name-filter <pat>    Case-insensitive filename filter (env: VCF_NAME_FILTER, default: *chr*)
  --bcftools-module <m>  Module to load (env: BCFTOOLS_MODULE, default: bcftools/1.18-gcc-12.3.0)
  --samples <file>       Optional sample list file for bcftools stats -s (env: SAMPLES_FILE)
  --dry-run              Print intended actions and exit 0 (env: DRY_RUN)
  --help, -h             Show this message and exit
EOF_USAGE
}

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

require_arg() {
    local opt="$1"
    local val="$2"
    if [ -z "${val}" ]; then
        echo "Missing value for ${opt}" >&2
        exit 1
    fi
}

wrap_name_filter() {
    local pattern="$1"
    if [[ "${pattern}" != *"*"* && "${pattern}" != *"?"* && "${pattern}" != *"["* ]]; then
        echo "*${pattern}*"
    else
        echo "${pattern}"
    fi
}

VCF_DIR="${VCF_DIR:-$PWD}"
WORK_DIR="${WORK_DIR:-}"
OUT_FILE="${OUT_FILE:-}"
VCF_NAME_FILTER="${VCF_NAME_FILTER:-*chr*}"
BCFTOOLS_MODULE="${BCFTOOLS_MODULE:-bcftools/1.18-gcc-12.3.0}"
SAMPLES_FILE="${SAMPLES_FILE:-}"
DRY_RUN="$(normalize_bool "${DRY_RUN:-false}")"

while [ "$#" -gt 0 ]; do
    case "$1" in
        --vcf-dir)
            require_arg "$1" "${2:-}"
            VCF_DIR="$2"
            shift 2
            ;;
        --work-dir)
            require_arg "$1" "${2:-}"
            WORK_DIR="$2"
            shift 2
            ;;
        --out)
            require_arg "$1" "${2:-}"
            OUT_FILE="$2"
            shift 2
            ;;
        --name-filter)
            require_arg "$1" "${2:-}"
            VCF_NAME_FILTER="$2"
            shift 2
            ;;
        --bcftools-module)
            require_arg "$1" "${2:-}"
            BCFTOOLS_MODULE="$2"
            shift 2
            ;;
        --samples)
            require_arg "$1" "${2:-}"
            SAMPLES_FILE="$2"
            shift 2
            ;;
        --dry-run)
            DRY_RUN="true"
            shift
            ;;
        --help|-h)
            usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1" >&2
            usage
            exit 1
            ;;
    esac
done

if ! VCF_DIR="$(cd "${VCF_DIR}" && pwd)"; then
    log_error "VCF_DIR does not exist or is not accessible: ${VCF_DIR}"
    exit 1
fi

if [ -z "${WORK_DIR}" ]; then
    WORK_DIR="${VCF_DIR}"
fi
mkdir -p "${WORK_DIR}"
WORK_DIR="$(cd "${WORK_DIR}" && pwd)"

if [ -z "${OUT_FILE}" ]; then
    OUT_FILE="${WORK_DIR}/snps_per_sample_counts.csv.gz"
fi
OUT_DIR="$(dirname "${OUT_FILE}")"
mkdir -p "${OUT_DIR}"
OUT_FILE="$(cd "${OUT_DIR}" && pwd)/$(basename "${OUT_FILE}")"

VCF_NAME_FILTER="$(wrap_name_filter "${VCF_NAME_FILTER}")"

mapfile -t VCF_FILES < <(find "${VCF_DIR}" -maxdepth 1 -type f -name "*.vcf.gz" -iname "${VCF_NAME_FILTER}" | LC_ALL=C sort)

if [ "${#VCF_FILES[@]}" -eq 0 ]; then
    log_error "No VCFs found in ${VCF_DIR} matching *.vcf.gz and filter ${VCF_NAME_FILTER}"
    exit 1
fi

log_info "Detected ${#VCF_FILES[@]} VCFs in ${VCF_DIR}"
log_info "Filters: bcftools view -m2 -M2 -v snps -f PASS | bcftools stats -s <samples>"

if [ "${DRY_RUN}" = "true" ]; then
    log_info "Dry-run enabled; no output will be created."
    log_info "VCF_DIR: ${VCF_DIR}"
    log_info "WORK_DIR: ${WORK_DIR}"
    log_info "OUT_FILE: ${OUT_FILE}"
    log_info "VCF_NAME_FILTER: ${VCF_NAME_FILTER}"
    log_info "BCFTOOLS_MODULE: ${BCFTOOLS_MODULE}"
    if [ -n "${SAMPLES_FILE}" ]; then
        log_info "SAMPLES_FILE: ${SAMPLES_FILE}"
    else
        log_info "SAMPLES_FILE: (all samples)"
    fi
    for vcf in "${VCF_FILES[@]}"; do
        log_info "Would process: ${vcf}"
    done
    exit 0
fi

if command -v module >/dev/null 2>&1; then
    module purge
    if ! module load "${BCFTOOLS_MODULE}" >/dev/null 2>&1; then
        log_error "Failed to load module: ${BCFTOOLS_MODULE}"
        exit 1
    fi
    log_info "Loaded module: ${BCFTOOLS_MODULE}"
else
    log_warn "module command not found; assuming bcftools is already in PATH"
fi

if ! command -v bcftools >/dev/null 2>&1; then
    log_error "bcftools not found in PATH after module load"
    exit 1
fi

SAMPLES_ARG="-"
if [ -n "${SAMPLES_FILE}" ]; then
    if [ ! -f "${SAMPLES_FILE}" ]; then
        log_error "Samples file does not exist: ${SAMPLES_FILE}"
        exit 1
    fi
    SAMPLES_FILE="$(cd "$(dirname "${SAMPLES_FILE}")" && pwd)/$(basename "${SAMPLES_FILE}")"
    SAMPLES_ARG="${SAMPLES_FILE}"
fi

TMP_DIR="$(mktemp -d)"
cleanup() {
    rm -rf "${TMP_DIR}"
}
trap cleanup EXIT

CHR_ROWS_TMP="${TMP_DIR}/chr_rows.csv"
TOTAL_ROWS_TMP="${TMP_DIR}/total_rows.csv"
: > "${CHR_ROWS_TMP}"
: > "${TOTAL_ROWS_TMP}"

FIRST_VCF="${VCF_FILES[0]}"
FIRST_SAMPLES="${TMP_DIR}/samples_first.txt"
if ! bcftools query -l "${FIRST_VCF}" > "${FIRST_SAMPLES}"; then
    log_error "Failed to read sample list from ${FIRST_VCF}"
    exit 1
fi

for vcf in "${VCF_FILES[@]:1}"; do
    if ! bcftools query -l "${vcf}" > "${TMP_DIR}/samples_this.txt"; then
        log_error "Failed to read sample list from ${vcf}"
        exit 1
    fi
    if ! cmp -s "${FIRST_SAMPLES}" "${TMP_DIR}/samples_this.txt"; then
        log_error "Sample list mismatch detected in ${vcf} (compared to ${FIRST_VCF})"
        exit 1
    fi
    rm -f "${TMP_DIR}/samples_this.txt"
done

SAMPLE_COUNT="$(wc -l < "${FIRST_SAMPLES}" | tr -d ' ')"
log_info "Sample count: ${SAMPLE_COUNT}"

for vcf in "${VCF_FILES[@]}"; do
    chr_label="$(basename "${vcf}")"
    chr_label="${chr_label%.vcf.gz}"
    log_info "Processing ${chr_label}"

    if ! bcftools view -m2 -M2 -v snps -f PASS -Ou "${vcf}" \
        | bcftools stats -s "${SAMPLES_ARG}" - \
        | awk -v chr_label="${chr_label}" -v source_vcf="${vcf}" -v OFS="," '
            BEGIN {
                header_seen = 0
                miss_idx = ""
            }
            $0 ~ /^# PSC/ {
                if ($0 !~ /\[[0-9]+\]/) {
                    next
                }
                header_seen = 1
                line = $0
                while (match(line, /\[([0-9]+)\]([A-Za-z0-9_]+)/, m)) {
                    idx[m[2]] = m[1] + 0
                    line = substr(line, RSTART + RLENGTH)
                }
                if (!("sample" in idx) || !("nRefHom" in idx) || !("nNonRefHom" in idx) || !("nHets" in idx)) {
                    print "ERROR: Missing required PSC columns in stats header for " source_vcf > "/dev/stderr"
                    exit 1
                }
                if ("nMissing" in idx) {
                    miss_idx = idx["nMissing"]
                } else if ("nMiss" in idx) {
                    miss_idx = idx["nMiss"]
                } else if ("nMissings" in idx) {
                    miss_idx = idx["nMissings"]
                } else {
                    print "ERROR: Missing nMissing PSC column in stats header for " source_vcf > "/dev/stderr"
                    exit 1
                }
                next
            }
            $1 == "PSC" {
                if (!header_seen) {
                    print "ERROR: PSC header not found in stats output for " source_vcf > "/dev/stderr"
                    exit 1
                }
                sample = $(idx["sample"])
                n_refhom = $(idx["nRefHom"]) + 0
                n_homalt = $(idx["nNonRefHom"]) + 0
                n_het = $(idx["nHets"]) + 0
                n_missing = $(miss_idx) + 0
                n_called = n_refhom + n_homalt + n_het
                n_total = n_called + n_missing
                print "CHR", chr_label, sample, n_called, n_missing, n_total, n_refhom, n_homalt, n_het, source_vcf
            }
            END {
                if (!header_seen) {
                    print "ERROR: PSC header not found in stats output for " source_vcf > "/dev/stderr"
                    exit 1
                }
            }
        ' >> "${CHR_ROWS_TMP}"; then
        log_error "Failed to process ${vcf}"
        exit 1
    fi
done

awk -F ',' -v OFS=',' '{
        sample = $3
        n_called[sample] += $4
        n_missing[sample] += $5
        n_total[sample] += $6
        n_refhom[sample] += $7
        n_homalt[sample] += $8
        n_het[sample] += $9
    }
    END {
        for (s in n_called) {
            print "TOTAL", "TOTAL", s, n_called[s], n_missing[s], n_total[s], n_refhom[s], n_homalt[s], n_het[s], ""
        }
    }
' "${CHR_ROWS_TMP}" > "${TOTAL_ROWS_TMP}"

{
    echo "scope,chr_label,sample,n_called,n_missing,n_total,n_refhom,n_homalt,n_het,source_vcf"
    LC_ALL=C sort -t ',' -k2,2 -k3,3 "${CHR_ROWS_TMP}"
    LC_ALL=C sort -t ',' -k3,3 "${TOTAL_ROWS_TMP}"
} | gzip -c > "${OUT_FILE}"

log_info "Wrote output to ${OUT_FILE}"
