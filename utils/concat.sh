#!/bin/bash --login
# =============================================================================
# Concatenate chromosome 1-17 VCFs on node-local storage.
#
# Usage:
#   sbatch utils/concat_chr01_chr17_noCommon_vcf.sh <vcf_dir> [output.vcf.gz]
#
# Example:
#   sbatch utils/concat_chr01_chr17_noCommon_vcf.sh \
#     /scratch/project/bigdata_apple/<dataset>/7.Consolidated_VCF \
#     /scratch/project/bigdata_apple/<dataset>/7.Consolidated_VCF/Chr01-Chr17_consolidated_noCommon.vcf.gz
#
# Input files are detected from non-recursive *.vcf.gz files whose basenames
# contain Chr1/Chr01 through Chr17 chromosome tokens.
#
# The output is created and indexed under $TMPDIR, then copied to the requested
# durable path only after bcftools completes successfully. $TMPDIR must have
# room for all 17 input files and the concatenated output.
# =============================================================================
#SBATCH --job-name=concat_noCommon_vcf
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1-00:00:00
#SBATCH --partition=general
#SBATCH --account=a_qaafi_cas
#SBATCH --output=slurm.concat_noCommon_vcf.%j.out
#SBATCH --error=slurm.concat_noCommon_vcf.%j.err
# =============================================================================

set -euo pipefail

usage() {
    cat >&2 <<'EOF_USAGE'
Usage: concat.sh <vcf_dir> [output.vcf.gz]

Arguments:
  vcf_dir        Directory containing chromosome 1-17 *.vcf.gz files.
  output.vcf.gz  Final compressed output VCF path.
                 Default: <vcf_dir>/Chr01-Chr17_consolidated_noCommon.vcf.gz

Environment:
  TMPDIR           Required node-local job temporary directory.
  BCFTOOLS_BIN     bcftools executable name or path (default: bcftools).
  BCFTOOLS_MODULE  Module loaded on Bunya when BCFTOOLS_BIN=bcftools
                   (default: bcftools/1.18-gcc-12.3.0).
EOF_USAGE
}

log_info() {
    printf '[INFO] %s\n' "$*"
}

log_error() {
    printf '[ERROR] %s\n' "$*" >&2
}

find_chromosome_vcf() {
    local chr_num="$1"
    local matches=()
    local candidate basename match_count

    shopt -s nullglob
    for candidate in "${VCF_DIR}"/*.vcf.gz; do
        [ -f "${candidate}" ] || continue
        [ "${candidate}" != "${OUTPUT_VCF:-}" ] || continue
        basename="$(basename "${candidate}")"
        # Skip range-named merged outputs from earlier concat runs.
        if [[ "${basename}" =~ [Cc][Hh][Rr]0*1[^[:digit:]][Cc][Hh][Rr]0*17([^[:digit:]]|$) ]]; then
            continue
        fi
        if [[ "${basename}" =~ (^|[^[:alnum:]])[Cc][Hh][Rr]0*${chr_num}([^[:digit:]]|$) ]]; then
            matches+=("${candidate}")
        fi
    done
    shopt -u nullglob

    match_count="${#matches[@]}"
    if [ "${match_count}" -eq 0 ]; then
        log_error "Required chromosome ${chr_num} VCF is missing in ${VCF_DIR}"
        exit 1
    fi
    if [ "${match_count}" -gt 1 ]; then
        log_error "Multiple VCFs matched chromosome ${chr_num}; expected exactly one:"
        printf '  %s\n' "${matches[@]}" >&2
        exit 1
    fi

    printf '%s\n' "${matches[0]}"
}

if [ "$#" -lt 1 ] || [ "$#" -gt 2 ]; then
    usage
    exit 1
fi

VCF_DIR="$1"
if ! VCF_DIR="$(cd "${VCF_DIR}" && pwd)"; then
    log_error "Input VCF directory does not exist or is not accessible: $1"
    exit 1
fi

OUTPUT_VCF="${2:-${VCF_DIR}/Chr01-Chr17_consolidated_noCommon.vcf.gz}"
if [[ "${OUTPUT_VCF}" != *.vcf.gz ]]; then
    log_error "Output must end in .vcf.gz: ${OUTPUT_VCF}"
    exit 1
fi
OUTPUT_DIR="$(dirname "${OUTPUT_VCF}")"
mkdir -p "${OUTPUT_DIR}"
OUTPUT_DIR="$(cd "${OUTPUT_DIR}" && pwd)"
OUTPUT_VCF="${OUTPUT_DIR}/$(basename "${OUTPUT_VCF}")"

if [ -z "${TMPDIR:-}" ] || [ ! -d "${TMPDIR}" ] || [ ! -w "${TMPDIR}" ]; then
    log_error "TMPDIR must be set to a writable node-local directory. Submit this script with sbatch."
    exit 1
fi

BCFTOOLS_BIN="${BCFTOOLS_BIN:-bcftools}"
BCFTOOLS_MODULE="${BCFTOOLS_MODULE:-bcftools/1.18-gcc-12.3.0}"
if [ "${BCFTOOLS_BIN}" = "bcftools" ] && command -v module >/dev/null 2>&1; then
    module purge
    module load "${BCFTOOLS_MODULE}"
fi
if ! command -v "${BCFTOOLS_BIN}" >/dev/null 2>&1; then
    log_error "bcftools executable not found: ${BCFTOOLS_BIN}"
    exit 1
fi

THREADS="${SLURM_CPUS_PER_TASK:-1}"
JOB_TAG="${SLURM_JOB_ID:-local}"
WORK_DIR="$(mktemp -d "${TMPDIR%/}/concat_noCommon_vcf.${JOB_TAG}.XXXXXX")"
cleanup() {
    rm -rf "${WORK_DIR}"
}
trap cleanup EXIT

INPUT_LIST="${WORK_DIR}/inputs.list"
: > "${INPUT_LIST}"

log_info "Detecting and staging chromosome 1-17 *.vcf.gz inputs to ${WORK_DIR}"
for chr_num in {1..17}; do
    source_vcf="$(find_chromosome_vcf "${chr_num}")"
    filename="$(basename "${source_vcf}")"
    if [ "${source_vcf}" = "${OUTPUT_VCF}" ]; then
        log_error "Output path would overwrite an input VCF: ${OUTPUT_VCF}"
        exit 1
    fi
    cp -p "${source_vcf}" "${WORK_DIR}/${filename}"
    printf '%s\n' "${WORK_DIR}/${filename}" >> "${INPUT_LIST}"
done

LOCAL_OUTPUT="${WORK_DIR}/$(basename "${OUTPUT_VCF}")"
log_info "Concatenating ordered inputs on TMPDIR into ${LOCAL_OUTPUT}"
"${BCFTOOLS_BIN}" concat \
    --threads "${THREADS}" \
    --file-list "${INPUT_LIST}" \
    --output-type z \
    --output "${LOCAL_OUTPUT}"

log_info "Creating CSI index on TMPDIR"
"${BCFTOOLS_BIN}" index --threads "${THREADS}" --csi --force "${LOCAL_OUTPUT}"

log_info "Copying completed VCF and CSI index to ${OUTPUT_DIR}"
cp -p "${LOCAL_OUTPUT}" "${OUTPUT_VCF}"
cp -p "${LOCAL_OUTPUT}.csi" "${OUTPUT_VCF}.csi"

log_info "Complete: ${OUTPUT_VCF}"
log_info "Index: ${OUTPUT_VCF}.csi"
