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

Repairs malformed VCF rows by ensuring at least 9 columns (CHROM..FORMAT) and
one genotype column per sample, filling missing values with ./.

Outputs:
  - Writes patched VCF (bgzip-compressed) to the given output path, or
    <input>.fixed.vcf.gz if not specified.
  - Prints a summary of total rows and how many were patched.
EOF
}

if [ $# -lt 1 ]; then
    usage
    exit 1
fi

INPUT="$1"
OUTPUT="${2:-${INPUT%.vcf.gz}.fixed.vcf.gz}"

# Ensure bcftools is available; load specific module if needed
BCFTOOLS_BIN="${BCFTOOLS_BIN:-bcftools}"
if ! command -v "${BCFTOOLS_BIN}" >/dev/null 2>&1; then
    if command -v module >/dev/null 2>&1; then
        module load bcftools/1.18-gcc-12.3.0 >/dev/null 2>&1 || true
    fi
fi
if ! command -v "${BCFTOOLS_BIN}" >/dev/null 2>&1; then
    echo "[FATAL] bcftools not found in PATH (after attempting module load bcftools/1.18-gcc-12.3.0)." >&2
    exit 1
fi

for tool in bgzip tabix; do
    if ! command -v "${tool}" >/dev/null 2>&1; then
        echo "[FATAL] ${tool} not found in PATH." >&2
        exit 1
    fi
done

if [ ! -f "${INPUT}" ]; then
    echo "[FATAL] Input VCF not found: ${INPUT}" >&2
    exit 1
fi

# Extract header to determine sample count
header_line=$(zcat "${INPUT}" 2>/dev/null | awk 'BEGIN{FS="\t"} /^#CHROM/ {print; exit}')
if [ -z "${header_line}" ]; then
    echo "[FATAL] #CHROM header not found in ${INPUT}" >&2
    exit 1
fi
header_fields=$(printf "%s" "${header_line}" | awk '{print NF}')
if [ "${header_fields}" -lt 9 ]; then
    echo "[FATAL] Header has fewer than 9 columns: ${header_fields}" >&2
    exit 1
fi
sample_count=$((header_fields - 9))
expected_cols=$((9 + sample_count))

total_rows=$(zcat "${INPUT}" 2>/dev/null | awk '!/^#/ {c++} END{print c+0}')
malformed_rows=$(zcat "${INPUT}" 2>/dev/null | awk -v expected="${expected_cols}" 'BEGIN{FS="\t"} !/^#/ && NF<expected {c++} END{print c+0}')

echo "[INFO] Input: ${INPUT}"
echo "[INFO] Samples: ${sample_count} | Expected columns: ${expected_cols}"
echo "[INFO] Total data rows: ${total_rows}"
echo "[INFO] Malformed rows (NF < ${expected_cols}): ${malformed_rows}"

tmp_base="${TMPDIR:-/tmp}"
tmp_out="$(mktemp "${tmp_base%/}/fixvcfXXXXXX")"
trap 'rm -f "${tmp_out}"' EXIT

zcat "${INPUT}" 2>/dev/null | \
awk -v expected="${expected_cols}" 'BEGIN{FS=OFS="\t"}
    /^#/ {print; next}
    {
        # Ensure at least CHROM..INFO (8 cols)
        if (NF < 8) {
            for (i = NF + 1; i <= 8; i++) {$i="."}
        }
        # Ensure FORMAT
        if (NF < 9) {
            $9 = "GT"
        }
        # Pad sample columns with ./.
        while (NF < expected) {
            $(NF + 1) = "./."
        }
        print
    }' > "${tmp_out}"

bgzip -c "${tmp_out}" > "${OUTPUT}"
tabix -f "${OUTPUT}"

# Verify patched file
if ! "${BCFTOOLS_BIN}" view -Ov -o /dev/null "${OUTPUT}" >/dev/null 2>&1; then
    echo "[FATAL] Patched VCF failed validation: ${OUTPUT}" >&2
    exit 1
fi

echo "[INFO] Patched VCF written: ${OUTPUT}"
patched_rows=${malformed_rows}
if [ "${total_rows}" -gt 0 ]; then
    percent=$(python - <<PY
total=${total_rows}
patched=${patched_rows}
print(f"{(patched/total*100):.2f}")
PY
)
else
    percent="0.00"
fi
echo "[INFO] Patched rows: ${patched_rows}/${total_rows} (${percent}%)"
