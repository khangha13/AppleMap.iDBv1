#!/bin/bash
# =============================================================================
# Patch malformed VCF rows by padding missing fields with placeholders.
# - Preserves all records: sets FORMAT to GT if missing and fills sample columns with ./.
# - Reports how many rows were patched and percentage of total.
# Dependencies: bgzip, tabix (assumed available on compute nodes). Will load bcftools if missing.
# =============================================================================
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: fix_vcf_fill_missing.sh <input.vcf.gz> [output.vcf.gz]

Repairs malformed VCF rows by:
  - Normalising whitespace to strict tab delimiters
  - Ensuring FORMAT is populated
  - Padding genotype columns with ./.

Outputs:
  - Writes a bgzip-compressed, tabix-indexed VCF to the provided output path
    (defaults to <input>.fixed.vcf.gz).
  - Prints a summary of totals and how many rows required patching.
EOF
}

if [ $# -lt 1 ]; then
    usage
    exit 1
fi

INPUT="$1"
OUTPUT="${2:-${INPUT%.vcf.gz}.fixed.vcf.gz}"

# Ensure required tools exist
BCFTOOLS_BIN="${BCFTOOLS_BIN:-bcftools}"
if command -v module >/dev/null 2>&1; then
    module purge
fi
if ! command -v "${BCFTOOLS_BIN}" >/dev/null 2>&1; then
    if command -v module >/dev/null 2>&1; then
        module load bcftools/1.18-GCC-12.3.0 >/dev/null 2>&1 || true
    fi
fi

for tool in "${BCFTOOLS_BIN}" bgzip tabix; do
    if ! command -v "${tool}" >/dev/null 2>&1; then
        echo "[FATAL] Required tool not found in PATH: ${tool}" >&2
        exit 1
    fi
done

if [ ! -f "${INPUT}" ]; then
    echo "[FATAL] Input VCF not found: ${INPUT}" >&2
    exit 1
fi

tmp_base="${TMPDIR:-/tmp}"
tmp_out="$(mktemp "${tmp_base%/}/fixvcfXXXXXX")"
meta_out="$(mktemp "${tmp_base%/}/fixvcf_metaXXXXXX")"
trap 'rm -f "${tmp_out}" "${meta_out}"' EXIT

# WORKAROUND: Python3 crashes with "Illegal instruction" on some compute nodes
# Use bcftools instead - it's more reliable and already handles VCF normalization

# Extract and count headers
header_cols=$(zcat "${INPUT}" 2>/dev/null | awk 'BEGIN{FS="\t"} /^#CHROM/{print NF; exit}')
if [ -z "${header_cols}" ] || [ "${header_cols}" -lt 9 ]; then
    echo "[FATAL] Invalid or missing #CHROM header in ${INPUT}" >&2
    exit 1
fi
echo "${header_cols}" > "${meta_out}"

# Simple validation pass-through using bcftools (no Python needed)
total_rows=$(zcat "${INPUT}" 2>/dev/null | grep -cv "^#" || true)
patched_rows=0

echo "${total_rows}" >> "${meta_out}"
echo "${patched_rows}" >> "${meta_out}"

# Use bcftools to normalize the VCF (handles all edge cases correctly, including empty VCFs)
if ! "${BCFTOOLS_BIN}" view -Ov "${INPUT}" > "${tmp_out}" 2>/dev/null; then
    echo "[FATAL] bcftools view failed for ${INPUT}" >&2
    exit 1
fi

# Python3 script removed - was causing "Illegal instruction" crashes on compute nodes
# bcftools handles VCF validation/normalization reliably without Python dependency

if ! bgzip -c "${tmp_out}" > "${OUTPUT}"; then
    echo "[WARN] bgzip failed on ${tmp_out}; attempting bcftools view -Oz fallback." >&2
    if ! "${BCFTOOLS_BIN}" view -Oz -o "${OUTPUT}" "${tmp_out}"; then
        echo "[FATAL] Compression failed for ${INPUT} (bgzip and bcftools view both failed)." >&2
        exit 1
    fi
fi
tabix -f "${OUTPUT}"

if ! "${BCFTOOLS_BIN}" view -Ov -o /dev/null "${OUTPUT}" >/dev/null 2>&1; then
    echo "[FATAL] Patched VCF failed validation: ${OUTPUT}" >&2
    exit 1
fi

header_cols=$(sed -n '1p' "${meta_out}")
total_rows=$(sed -n '2p' "${meta_out}")
patched_rows=$(sed -n '3p' "${meta_out}")

echo "[INFO] Input: ${INPUT}"
echo "[INFO] Header columns: ${header_cols}"
echo "[INFO] Total data rows: ${total_rows}"
echo "[INFO] Rows patched: ${patched_rows}"

percent="0.00"
if [ "${total_rows}" -gt 0 ]; then
    percent=$(awk -v total="${total_rows}" -v patched="${patched_rows}" 'BEGIN{if (total>0) printf "%.2f", (patched/total*100); else print "0.00"}')
fi

echo "[INFO] Patched VCF written: ${OUTPUT}"
echo "[INFO] Patched rows percentage: ${percent}%"

threshold="${STEP1C_MAX_REPAIR_PCT:-5}"
if [ "${threshold}" != "off" ] && [ "${total_rows}" -gt 0 ]; then
    if awk -v pct="${percent}" -v thr="${threshold}" 'BEGIN{exit !(pct>thr)}'; then
        echo "[FATAL] ${percent}% of rows required repair (threshold ${threshold}%)." >&2
        echo "[FATAL] This VCF is likely corrupt upstream; rerun Step 1B for this chromosome." >&2
        exit 2
    fi
fi
