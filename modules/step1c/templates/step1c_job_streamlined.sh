#!/bin/bash
# STEP 1C IMPUTATION JOB TEMPLATE (Streamlined Version)
# Removes redundant validation loops now that bgzip incompatibility is the known root cause
set -euo pipefail

DATASET_PATH="$1"
DATASET_NAME="$(basename "${DATASET_PATH}")"
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

STEP1C_RUN_ID="$(date +%Y%m%d_%H%M%S)"
if command -v resolve_log_root >/dev/null 2>&1; then
    STEP1C_DEBUG_DIR="$(resolve_log_root "${DATASET_NAME}" "artifacts/step1c_debug/${STEP1C_RUN_ID}")"
else
    STEP1C_DEBUG_DIR="${LOG_BASE_PATH%/}/${DATASET_NAME}/artifacts/step1c_debug/${STEP1C_RUN_ID}"
fi
STEP1C_PERSIST_ARTIFACTS="${STEP1C_PERSIST_ARTIFACTS:-true}"
STEP1C_COPY_BEAGLE_VCFS="${STEP1C_DEBUG_BEAGLE:-false}"
STEP1C_ARTIFACTS_NOTED="false"

persist_step1c_artifacts() {
    if [ "${STEP1C_PERSIST_ARTIFACTS}" != "true" ] || [ -z "${WORK_TMPDIR:-}" ]; then
        return
    fi
    mkdir -p "${STEP1C_DEBUG_DIR}"
    local patterns=("*.log" "*.txt")
    local found=false
    for pattern in "${patterns[@]}"; do
        while IFS= read -r artifact; do
            [ -z "${artifact}" ] && continue
            rsync -rt "${artifact}" "${STEP1C_DEBUG_DIR}/" >/dev/null 2>&1 || true
            found=true
        done < <(find "${WORK_TMPDIR}" -maxdepth 1 -type f -name "${pattern}" -print 2>/dev/null)
    done
    if [ "${STEP1C_COPY_BEAGLE_VCFS}" = "true" ]; then
        find "${WORK_TMPDIR}" -maxdepth 1 -type f \( -name "*.vcf.gz" -o -name "*.vcf.gz.tbi" -o -name "*.vcf.gz.csi" \) -exec rsync -rt {} "${STEP1C_DEBUG_DIR}/" \; >/dev/null 2>&1 || true
    fi
    if [ "${STEP1C_ARTIFACTS_NOTED}" = "false" ] && { [ "${found}" = "true" ] || [ "${STEP1C_COPY_BEAGLE_VCFS}" = "true" ]; }; then
        printf '%s\tPersisted artifacts from %s\n' "$(date +%Y-%m-%dT%H:%M:%S)" "${WORK_TMPDIR}" >> "${STEP1C_DEBUG_DIR}/ARTIFACTS.log"
        STEP1C_ARTIFACTS_NOTED="true"
    fi
}
trap 'persist_step1c_artifacts' EXIT

log_info "Using working directory: ${WORK_TMPDIR}"
log_info "Step 1C artifacts will be copied to ${STEP1C_DEBUG_DIR}"

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
LOCAL_REF_FAI="${LOCAL_REF}.fai"
numeric_contig_header="${WORK_TMPDIR}/contigs_numeric.header.txt"

VCF_PREFIXES=()
BEAGLE_VCF_PATHS=()

# Process each input VCF: filter, normalize, rename contigs, sort
# Removed redundant auto-fix and validation loops (root cause was bgzip, not VCF format)
while IFS= read -r vcf_path; do
    [ -z "${vcf_path}" ] && continue
    [ -f "${vcf_path}" ] || error_exit "Input VCF not found: ${vcf_path}"
    
    base_name="$(basename "${vcf_path}")"
    log_info "Processing VCF: ${vcf_path}"
    
    # Step 1: Filter to biallelic SNPs and normalize against reference
    filtered_base="${base_name%.vcf.gz}.filtered.vcf.gz"
    filtered_path="${WORK_TMPDIR}/${filtered_base}"
    log_info "Filtering/normalizing SNPs for ${base_name}"
    if ! "${BCFTOOLS_BIN}" norm -m -any -f "${LOCAL_REF}" "${vcf_path}" -Ou | \
         "${BCFTOOLS_BIN}" view -v snps -m2 -M2 -Oz -o "${filtered_path}"; then
        error_exit "bcftools filter/normalize failed for ${base_name}"
    fi
    if ! "${BCFTOOLS_BIN}" index -f -c "${filtered_path}"; then
        log_warn "Failed to index ${filtered_base}"
    fi

    # Step 2: Rename contigs (Chr01 -> 1, Chr02 -> 2, etc.)
    renamed_base="${filtered_base%.vcf.gz}.renamed.vcf.gz"
    renamed_path="${WORK_TMPDIR}/${renamed_base}"
    log_info "Renaming contigs (Chr01->1) in ${filtered_base}"
    if ! "${BCFTOOLS_BIN}" annotate --rename-chrs "${mapping_file}" -Oz -o "${renamed_path}" "${filtered_path}"; then
        error_exit "bcftools annotate (rename) failed for ${filtered_base}"
    fi
    if ! "${BCFTOOLS_BIN}" index -f -c "${renamed_path}"; then
        log_warn "Failed to index ${renamed_base}"
    fi

    # Step 3: Inject numeric contig headers if reference FAI exists
    final_path="${renamed_path}"
    final_base="${renamed_base}"
    if [ -f "${LOCAL_REF_FAI}" ]; then
        if [ ! -s "${numeric_contig_header}" ]; then
            awk 'NR==FNR {map[$1]=$2; next} {if ($1 in map) printf "##contig=<ID=%s,length=%s>\n", map[$1], $2}' \
                "${mapping_file}" "${LOCAL_REF_FAI}" > "${numeric_contig_header}"
        fi
        if [ -s "${numeric_contig_header}" ]; then
            header_base="${renamed_base%.vcf.gz}.hdr.vcf.gz"
            header_path="${WORK_TMPDIR}/${header_base}"
            log_info "Injecting numeric contig header into ${renamed_base}"
            if ! "${BCFTOOLS_BIN}" annotate -h "${numeric_contig_header}" -Oz -o "${header_path}" "${renamed_path}"; then
                error_exit "bcftools annotate (add contig header) failed for ${renamed_base}"
            fi
            if ! "${BCFTOOLS_BIN}" index -f -c "${header_path}"; then
                log_warn "Failed to index ${header_base}"
            fi
            final_path="${header_path}"
            final_base="${header_base}"
        fi
    fi

    # Step 4: Sort VCF for Beagle
    sorted_base="${final_base%.vcf.gz}.sorted.beagle.vcf.gz"
    sorted_path="${WORK_TMPDIR}/${sorted_base}"
    log_info "Sorting VCF for Beagle: ${final_base}"
    if ! "${BCFTOOLS_BIN}" sort "${final_path}" -Oz -o "${sorted_path}"; then
        error_exit "bcftools sort failed for ${final_base}"
    fi
    if ! "${BCFTOOLS_BIN}" index -f -c "${sorted_path}"; then
        log_warn "Failed to index ${sorted_base}"
    fi

    # Quick sanity check: ensure VCF has proper structure
    if ! "${BCFTOOLS_BIN}" view -h "${sorted_path}" | grep -q "^#CHROM"; then
        error_exit "VCF missing #CHROM header line: ${sorted_base}"
    fi

    VCF_PREFIXES+=("${sorted_base}")
    BEAGLE_VCF_PATHS+=("${sorted_path}")
done < "${VCF_MANIFEST}"

LOCAL_MANIFEST="${WORK_TMPDIR}/vcf_manifest.txt"

# WORKAROUND: Beagle 5.4.22Jul22.46e has compatibility issues with certain bgzip formats
# Decompress VCFs to plain text for Beagle (verified working with uncompressed VCFs)
log_info "Decompressing VCFs for Beagle (bgzip format compatibility workaround)..."
UNCOMPRESSED_PREFIXES=()
for path in "${BEAGLE_VCF_PATHS[@]}"; do
    [ -f "${path}" ] || continue
    base_name="$(basename "${path}")"
    uncompressed_base="${base_name%.vcf.gz}.vcf"
    uncompressed_path="${WORK_TMPDIR}/${uncompressed_base}"
    log_info "Decompressing: ${base_name} -> ${uncompressed_base}"
    if ! zcat "${path}" > "${uncompressed_path}"; then
        error_exit "Failed to decompress ${path} for Beagle"
    fi
    UNCOMPRESSED_PREFIXES+=("${uncompressed_path}")
done

# Update manifest to point to uncompressed files
printf '%s\n' "${UNCOMPRESSED_PREFIXES[@]}" > "${LOCAL_MANIFEST}"
log_info "Prepared $(echo "${UNCOMPRESSED_PREFIXES[@]}" | wc -w) uncompressed VCF files for Beagle"

local_map_arg=""
if [ -n "${GENE_MAP_FILE}" ] && [ -f "${GENE_MAP_FILE}" ]; then
    rsync -rhivPt "${GENE_MAP_FILE}" "${WORK_TMPDIR}/"
    local_map_arg="map=${WORK_TMPDIR}/$(basename "${GENE_MAP_FILE}")"
fi

cd "${WORK_TMPDIR}"

JAVA_MEM="-Xmx${MEMORY_GB}g"

log_info "Launching Beagle imputation..."

# Minimum variants required for Beagle phasing (configurable)
STEP1C_MIN_VARIANTS="${STEP1C_MIN_VARIANTS:-10}"

# Process each chromosome separately (Beagle does not support manifest files in gt= parameter)
for uncompressed_vcf in "${UNCOMPRESSED_PREFIXES[@]}"; do
    [ -f "${uncompressed_vcf}" ] || continue
    
    # Extract chromosome name for output prefix
    chr_prefix="${uncompressed_vcf%.vcf}"
    chr_prefix="${chr_prefix##*/}"
    
    # Count variants (non-header lines)
    variant_count=$(grep -cv "^#" "${uncompressed_vcf}" || echo "0")
    
    if [ "${variant_count}" -lt "${STEP1C_MIN_VARIANTS}" ]; then
        log_warn "Skipping Beagle for ${uncompressed_vcf}: only ${variant_count} variants (min: ${STEP1C_MIN_VARIANTS})"
        log_info "Copying input as 'phased' output (insufficient variants for phasing statistics)"
        
        # Compress and copy as "phased" output
        bgzip -c "${uncompressed_vcf}" > "${chr_prefix}_phased.vcf.gz"
        if ! "${BCFTOOLS_BIN}" index -f -c "${chr_prefix}_phased.vcf.gz"; then
            log_warn "Failed to index ${chr_prefix}_phased.vcf.gz"
        fi
        
        # Create a minimal log explaining the skip
        {
            echo "Beagle skipped for ${uncompressed_vcf}"
            echo "Reason: Only ${variant_count} variants (minimum required: ${STEP1C_MIN_VARIANTS})"
            echo "Output: Input VCF copied as-is (no phasing applied)"
        } > "${chr_prefix}_phased.log"
        
        continue
    fi
    
    log_info "Running Beagle on ${uncompressed_vcf} (${variant_count} variants)..."
    
    if ! java ${JAVA_MEM} -jar "${BEAGLE_JAR}" \
        gt="${uncompressed_vcf}" \
        out="${chr_prefix}_phased" \
        nthreads="${THREADS}" \
        impute="${IMPUTE_FLAG}" \
        window=3 overlap=0.3 ne=100000 seed=2025 \
        ${local_map_arg}; then
        error_exit "Beagle failed for ${uncompressed_vcf}"
    fi
    
    log_info "Beagle completed for ${uncompressed_vcf}"
done

mkdir -p "${OUTPUT_DIR}"
rsync -rhivPt "${WORK_TMPDIR}"/*_phased.vcf.gz "${OUTPUT_DIR}/" || log_warn "No phased VCF outputs found"
rsync -rhivPt "${WORK_TMPDIR}"/*_phased.log "${OUTPUT_DIR}/" || log_warn "No Beagle log files found"

log_info "Beagle imputation completed. Results copied to ${OUTPUT_DIR}"

