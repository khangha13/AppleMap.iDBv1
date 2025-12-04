#!/bin/bash
# STEP 1C IMPUTATION JOB TEMPLATE
set -euo pipefail

DATASET_PATH="$1"
VCF_MANIFEST="$2"
OUTPUT_DIR="$3"
REFERENCE_FASTA="$4"
GENE_MAP_FILE="$5"
OUTPUT_PREFIX="$6"
THREADS="$7"
MEMORY_GB="$8"
IMPUTE_FLAG="${9:-false}"

# Decode sentinel for omitted gene map (preserves argument positions)
if [ "${GENE_MAP_FILE}" = "__NO_GENE_MAP__" ]; then
    GENE_MAP_FILE=""
fi

impute_flag_lc="$(echo "${IMPUTE_FLAG}" | tr '[:upper:]' '[:lower:]')"
case "${impute_flag_lc}" in
    true|1|yes|y)
        IMPUTE_FLAG="true"
        ;;
    false|0|no|n|"")
        IMPUTE_FLAG="false"
        ;;
    *)
        echo "[step1c_job] ⚠️  Unrecognized impute flag '${IMPUTE_FLAG}', defaulting to false (phasing-only)." >&2
        IMPUTE_FLAG="false"
        ;;
esac

if [ -n "${PIPELINE_ROOT:-}" ]; then
    if [ -d "${PIPELINE_ROOT}" ]; then
        PIPELINE_ROOT="$(cd "${PIPELINE_ROOT}" && pwd)"
    else
        echo "[step1c_job] ⚠️  Provided PIPELINE_ROOT (${PIPELINE_ROOT}) does not exist; falling back to template-relative path." >&2
        PIPELINE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
    fi
else
    PIPELINE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
fi
export PIPELINE_ROOT
STEP1C_MODULE_DIR="${PIPELINE_ROOT}/modules/step1c"
FIX_VCF_SCRIPT="${STEP1C_MODULE_DIR}/bin/fix_vcf_fill_missing.sh"

source "${PIPELINE_ROOT}/config/pipeline_config.sh"
source "${PIPELINE_ROOT}/lib/logging.sh"
source "${STEP1C_MODULE_DIR}/lib/functions.sh"

init_logging "step1c" "job"
log_info "Beagle impute flag: ${IMPUTE_FLAG} (false = phasing-only)"

# Load Beagle (fixed version for pipeline reproducibility)
module load beagle/5.4.22jul22.46e-java-11 >/dev/null 2>&1 || log_warn "Unable to load beagle/5.4.22jul22.46e-java-11 module; ensure beagle.jar is accessible."

BEAGLE_JAR="${BEAGLE_JAR:-${EBROOTBEAGLE:-}/beagle.jar}"
if [ -z "${BEAGLE_JAR}" ] || [ ! -f "${BEAGLE_JAR}" ]; then
    error_exit "Beagle jar not found. Set BEAGLE_JAR or load the beagle module."
fi

# Ensure bcftools is available (module load if needed)
BCFTOOLS_BIN="${BCFTOOLS_BIN:-bcftools}"
if ! command -v "${BCFTOOLS_BIN}" >/dev/null 2>&1; then
    module load bcftools/1.18-gcc-12.3.0 >/dev/null 2>&1 || log_warn "Unable to load bcftools/1.18-gcc-12.3.0; assuming bcftools is already available on compute node."
fi
if ! command -v "${BCFTOOLS_BIN}" >/dev/null 2>&1; then
    error_exit "bcftools not found in PATH. Install or load bcftools/1.18-gcc-12.3.0."
fi

WORK_TMPDIR="${TMPDIR:-$(mktemp -d "${SCRATCH_BASE_PATH%/}/step1c_XXXXXX")}"
mkdir -p "${WORK_TMPDIR}"

log_info "Using working directory: ${WORK_TMPDIR}"

mapping_file="${WORK_TMPDIR}/chr_rename_map.txt"
if [ ! -f "${mapping_file}" ]; then
    for i in $(seq -w 1 17); do
        printf "Chr%s\t%d\n" "${i}" "${i#0}" >> "${mapping_file}"
    done
fi

rsync -rhivPt "${REFERENCE_FASTA}" "${WORK_TMPDIR}/"
rsync -rhivPt "${REFERENCE_FASTA}.fai" "${WORK_TMPDIR}/" || true
rsync -rhivPt "${REFERENCE_FASTA%.*}.dict" "${WORK_TMPDIR}/" || true
LOCAL_REF="${WORK_TMPDIR}/$(basename "${REFERENCE_FASTA}")"

VCF_PREFIXES=()
while IFS= read -r vcf_path; do
    [ -z "${vcf_path}" ] && continue
    base_name="$(basename "${vcf_path}")"
    fixed_base="${base_name%.vcf.gz}.fixed.vcf.gz"
    fixed_path="${WORK_TMPDIR}/${fixed_base}"
    if [ ! -f "${FIX_VCF_SCRIPT}" ]; then
        error_exit "Auto-fix script not found: ${FIX_VCF_SCRIPT}"
    fi
    log_info "Auto-fixing VCF: ${vcf_path} -> ${fixed_path}"
    if ! bash "${FIX_VCF_SCRIPT}" "${vcf_path}" "${fixed_path}"; then
        error_exit "Auto-fix failed for ${vcf_path}"
    fi
    renamed_base="${fixed_base%.vcf.gz}.renamed.vcf.gz"
    renamed_path="${WORK_TMPDIR}/${renamed_base}"
    log_info "Renaming contigs (Chr01->1) in ${fixed_base}"
    if ! "${BCFTOOLS_BIN}" annotate --rename-chrs "${mapping_file}" -Oz -o "${renamed_path}" "${fixed_path}"; then
        error_exit "bcftools annotate (rename) failed for ${fixed_base}"
    fi
    tabix -f "${renamed_path}" || log_warn "Failed to index ${renamed_base}"

    # Filter to biallelic SNPs and normalize against reference
    filtered_base="${fixed_base%.vcf.gz}.filtered.vcf.gz"
    filtered_path="${WORK_TMPDIR}/${filtered_base}"
    log_info "Filtering/normalizing SNPs for ${renamed_base}"
    if ! "${BCFTOOLS_BIN}" norm -m -any -f "${LOCAL_REF}" "${renamed_path}" -Ou | \
         "${BCFTOOLS_BIN}" view -v snps -m2 -M2 -Oz -o "${filtered_path}"; then
        error_exit "bcftools filter/normalize failed for ${renamed_base}"
    fi
    tabix -f "${filtered_path}" || log_warn "Failed to index ${filtered_base}"

    local_name="${filtered_base}"
    VCF_PREFIXES+=("${local_name}")
done < "${VCF_MANIFEST}"

LOCAL_MANIFEST="${WORK_TMPDIR}/vcf_manifest.txt"
printf '%s\n' "${VCF_PREFIXES[@]}" > "${LOCAL_MANIFEST}"

local_map_arg=""
if [ -n "${GENE_MAP_FILE}" ] && [ -f "${GENE_MAP_FILE}" ]; then
    rsync -rhivPt "${GENE_MAP_FILE}" "${WORK_TMPDIR}/"
    local_map_arg="map=${WORK_TMPDIR}/$(basename "${GENE_MAP_FILE}")"
fi

cd "${WORK_TMPDIR}"

JAVA_MEM="-Xmx${MEMORY_GB}g"

log_info "Launching Beagle imputation..."

java ${JAVA_MEM} -jar "${BEAGLE_JAR}" \
    gt="${LOCAL_MANIFEST}" \
    out="${WORK_TMPDIR}/${OUTPUT_PREFIX}" \
    nthreads="${THREADS}" \
    impute="${IMPUTE_FLAG}" \
    window=3 overlap=0.3 ne=100000 seed=2025 \
    ${local_map_arg}

mkdir -p "${OUTPUT_DIR}"
rsync -rhivPt "${WORK_TMPDIR}/${OUTPUT_PREFIX}"* "${OUTPUT_DIR}/"

log_info "Beagle imputation completed. Results copied to ${OUTPUT_DIR}"
