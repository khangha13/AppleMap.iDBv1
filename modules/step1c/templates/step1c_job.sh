#!/bin/bash
# STEP 1C IMPUTATION JOB TEMPLATE
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
FIX_VCF_SCRIPT="${STEP1C_MODULE_DIR}/bin/fix_vcf_fill_missing.sh"
if [ -z "${STEP1C_MAX_REPAIR_PCT:-}" ]; then
    export STEP1C_MAX_REPAIR_PCT=5
else
    export STEP1C_MAX_REPAIR_PCT
fi

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
STEP1C_DEBUG_DIR="${LOG_BASE_PATH%/}/${DATASET_NAME}/step1c_debug/${STEP1C_RUN_ID}"
STEP1C_PERSIST_ARTIFACTS="${STEP1C_PERSIST_ARTIFACTS:-true}"
STEP1C_COPY_BEAGLE_VCFS="${STEP1C_DEBUG_BEAGLE:-false}"
STEP1C_ARTIFACTS_NOTED="false"

persist_step1c_artifacts() {
    if [ "${STEP1C_PERSIST_ARTIFACTS}" != "true" ] || [ -z "${WORK_TMPDIR:-}" ]; then
        return
    fi
    mkdir -p "${STEP1C_DEBUG_DIR}"
    local patterns=("*.validate.log" "*.dropped*.txt")
    local found=false
    for pattern in "${patterns[@]}"; do
        while IFS= read -r artifact; do
            [ -z "${artifact}" ] && continue
            rsync -rt "${artifact}" "${STEP1C_DEBUG_DIR}/" >/dev/null 2>&1 || true
            found=true
        done < <(find "${WORK_TMPDIR}" -maxdepth 1 -type f -name "${pattern}" -print 2>/dev/null)
    done
    if [ "${STEP1C_COPY_BEAGLE_VCFS}" = "true" ]; then
        find "${WORK_TMPDIR}" -maxdepth 1 -type f -name "*.beagle.vcf.gz*" -exec rsync -rt {} "${STEP1C_DEBUG_DIR}/" \; >/dev/null 2>&1 || true
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
    # Filter to biallelic SNPs and normalize against reference (using original contig names)
    filtered_base="${fixed_base%.vcf.gz}.filtered.vcf.gz"
    filtered_path="${WORK_TMPDIR}/${filtered_base}"
    log_info "Filtering/normalizing SNPs for ${fixed_base}"
    if ! "${BCFTOOLS_BIN}" norm -m -any -f "${LOCAL_REF}" "${fixed_path}" -Ou | \
         "${BCFTOOLS_BIN}" view -v snps -m2 -M2 -Oz -o "${filtered_path}"; then
        error_exit "bcftools filter/normalize failed for ${fixed_base}"
    fi
    tabix -f "${filtered_path}" || log_warn "Failed to index ${filtered_base}"

    # Rename contigs to numeric and inject corresponding header lines
    renamed_base="${filtered_base%.vcf.gz}.renamed.vcf.gz"
    renamed_path="${WORK_TMPDIR}/${renamed_base}"
    log_info "Renaming contigs (Chr01->1) in ${filtered_base}"
    if ! "${BCFTOOLS_BIN}" annotate --rename-chrs "${mapping_file}" -Oz -o "${renamed_path}" "${filtered_path}"; then
        error_exit "bcftools annotate (rename) failed for ${filtered_base}"
    fi
    tabix -f "${renamed_path}" || log_warn "Failed to index ${renamed_base}"

    header_path="${renamed_path}"
    header_base="${renamed_base}"
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
            tabix -f "${header_path}" || log_warn "Failed to index ${header_base}"
        fi
    fi

    # Final safety pass: re-pad/fix after header edits to catch any residual malformed rows
    final_base="${header_base%.vcf.gz}.final.vcf.gz"
    final_path="${WORK_TMPDIR}/${final_base}"
    log_info "Final padding/validation for ${header_base}"
    if ! bash "${FIX_VCF_SCRIPT}" "${header_path}" "${final_path}"; then
        error_exit "Final auto-fix failed for ${header_base}"
    fi

    # Drop any remaining rows with insufficient columns and log count
    clean_base="${final_base%.vcf.gz}.clean.vcf.gz"
    clean_path="${WORK_TMPDIR}/${clean_base}"
    dropped_log="${WORK_TMPDIR}/${final_base%.vcf.gz}.dropped.txt"
    validate_log="${WORK_TMPDIR}/${final_base%.vcf.gz}.validate.txt"
    expected_cols=$(zcat "${final_path}" 2>/dev/null | awk 'BEGIN{FS="\t"} /^#CHROM/{print NF; exit}')
    if [ -z "${expected_cols}" ]; then
        expected_cols=9
    fi
    log_info "Dropping rows with <${expected_cols} columns in ${final_base}"
    zcat "${final_path}" 2>/dev/null | awk -v expected="${expected_cols}" 'BEGIN{FS="\t"; OFS="\t"}
        /^#/ {print; next}
        {
            if (NF < expected) {
                print > "'"${dropped_log}"'"
                next
            }
            print
        }' | bgzip > "${clean_path}"
    tabix -f "${clean_path}" || log_warn "Failed to index ${clean_base}"
    if [ -s "${dropped_log}" ]; then
        log_warn "Dropped malformed rows from ${final_base}; details in ${dropped_log}"
    fi
    # Final validation: abort if any record is still short
    if ! zcat "${clean_path}" 2>/dev/null | awk -v expected="${expected_cols}" 'BEGIN{FS="\t"} /^#CHROM/{next} !/^#/ && NF<expected {print "Line " NR " has " NF " fields"; exit 1}'; then
        error_exit "Final VCF still contains rows with <${expected_cols} columns: ${clean_base}. See ${dropped_log} for dropped rows."
    fi

    # Sort the cleaned VCF by coordinate to match Beagle expectations
    sorted_base="${clean_base%.vcf.gz}.sorted.vcf.gz"
    sorted_path="${WORK_TMPDIR}/${sorted_base}"
    log_info "Sorting VCF for Beagle: ${clean_base}"
    if ! "${BCFTOOLS_BIN}" sort "${clean_path}" -Oz -o "${sorted_path}"; then
        error_exit "bcftools sort failed for ${clean_base}"
    fi
    tabix -f "${sorted_path}" || log_warn "Failed to index ${sorted_base}"

    # Strict validation: log and abort if any record has fewer columns than header
    expected_sorted=$(zcat "${sorted_path}" 2>/dev/null | awk 'BEGIN{FS="\t"} /^#CHROM/{print NF; exit}')
    if [ -z "${expected_sorted}" ]; then
        expected_sorted=9
    fi
    validate_report="${WORK_TMPDIR}/${sorted_base%.vcf.gz}.validate.txt"
    if ! zcat "${sorted_path}" 2>/dev/null | awk -v expected="${expected_sorted}" 'BEGIN{FS="\t"}
        /^#/ {next}
        NF<expected {print "Line " NR " has " NF " fields: "$0; exit 1}
    '; then
        {
            echo "Validation failed for ${sorted_base}"
            echo "Expected columns: ${expected_sorted}"
            zcat "${sorted_path}" 2>/dev/null | awk -v expected="${expected_sorted}" 'BEGIN{FS="\t"}
                /^#/ {next}
                NF<expected {print "Line " NR " has " NF " fields: "$0; exit}
            '
        } > "${validate_report}"
        error_exit "Final VCF still malformed: ${sorted_base}. See ${validate_report}"
    fi
    # Drop any residual short records and reindex (belt-and-braces)
    beagle_base="${sorted_base%.vcf.gz}.beagle.vcf.gz"
    beagle_path="${WORK_TMPDIR}/${beagle_base}"
    dropped_final_log="${WORK_TMPDIR}/${sorted_base%.vcf.gz}.dropped_final.txt"
    zcat "${sorted_path}" 2>/dev/null | awk -v expected="${expected_sorted}" 'BEGIN{FS="\t"; OFS="\t"}
        /^#/ {print; next}
        {
            if (NF < expected) {
                print > "'"${dropped_final_log}"'"
                next
            }
            print
        }' | bgzip > "${beagle_path}"
    tabix -f "${beagle_path}" || log_warn "Failed to index ${beagle_base}"
    if [ -s "${dropped_final_log}" ]; then
        log_warn "Dropped residual malformed rows from ${sorted_base}; details in ${dropped_final_log}"
    fi

    local_name="${beagle_base}"
    VCF_PREFIXES+=("${local_name}")
    BEAGLE_VCF_PATHS+=("${beagle_path}")
done < "${VCF_MANIFEST}"

LOCAL_MANIFEST="${WORK_TMPDIR}/vcf_manifest.txt"

VALIDATOR_SCRIPT="${STEP1C_MODULE_DIR}/bin/debug_beagle_vcf.sh"
if [ -x "${VALIDATOR_SCRIPT}" ]; then
    for path in "${BEAGLE_VCF_PATHS[@]}"; do
        [ -f "${path}" ] || continue
        validate_log="${path%.vcf.gz}.validate.log"
        log_info "Validating Beagle input: ${path}"
        if ! bash "${VALIDATOR_SCRIPT}" "${path}" > "${validate_log}" 2>&1; then
            error_exit "Beagle input validation failed: ${path}. See ${validate_log}"
        fi
    done
else
    log_warn "Beagle validator not found or not executable: ${VALIDATOR_SCRIPT}"
fi

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
    UNCOMPRESSED_PREFIXES+=("${uncompressed_base}")
done

# Update manifest to point to uncompressed files
printf '%s\n' "${UNCOMPRESSED_PREFIXES[@]}" > "${LOCAL_MANIFEST}"
log_info "Updated manifest with $(echo "${UNCOMPRESSED_PREFIXES[@]}" | wc -w) uncompressed VCF files"

local_map_arg=""
if [ -n "${GENE_MAP_FILE}" ] && [ -f "${GENE_MAP_FILE}" ]; then
    rsync -rhivPt "${GENE_MAP_FILE}" "${WORK_TMPDIR}/"
    local_map_arg="map=${WORK_TMPDIR}/$(basename "${GENE_MAP_FILE}")"
fi

cd "${WORK_TMPDIR}"

JAVA_MEM="-Xmx${MEMORY_GB}g"

log_info "Launching Beagle imputation..."

# Process each chromosome separately (Beagle does not support manifest files in gt= parameter)
for uncompressed_vcf in "${UNCOMPRESSED_PREFIXES[@]}"; do
    [ -f "${uncompressed_vcf}" ] || continue
    
    # Extract chromosome name for output prefix
    chr_prefix="${uncompressed_vcf%.vcf}"
    chr_prefix="${chr_prefix##*/}"  # Remove path
    
    log_info "Running Beagle on ${uncompressed_vcf}..."
    
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
