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
if ! command -v "${BCFTOOLS_BIN}" >/dev/null 2>&1; then
    if command -v module >/dev/null 2>&1; then
        module load bcftools/1.18-gcc-12.3.0 >/dev/null 2>&1 || true
    fi
fi

for tool in "${BCFTOOLS_BIN}" bgzip tabix python3; do
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

python3 - "$INPUT" "$tmp_out" "$meta_out" <<'PY'
import gzip
import sys

input_path, output_path, meta_path = sys.argv[1:4]

header_fields = None
total_rows = 0
patched_rows = 0

with gzip.open(input_path, "rt", encoding="utf-8", errors="replace") as fin, \
        open(output_path, "w", encoding="utf-8") as fout:
    for raw in fin:
        if raw.startswith("#"):
            line = raw.rstrip("\r\n")
            fout.write(line + "\n")
            if raw.startswith("#CHROM"):
                header_fields = line.split("\t")
            continue

        if header_fields is None:
            raise SystemExit("[FATAL] #CHROM header not found before data lines")

        total_rows += 1
        expected_cols = len(header_fields)
        tokens = raw.rstrip("\r\n")

        if "\t" in tokens:
            fields = tokens.split("\t")
        else:
            fields = tokens.split()

        if len(fields) != expected_cols:
            patched_rows += 1

        if len(fields) < expected_cols:
            fields.extend([""] * (expected_cols - len(fields)))
        elif len(fields) > expected_cols:
            fields = fields[:expected_cols]

        # Normalise critical columns
        fields[0] = fields[0] or "."
        fields[1] = fields[1] or "0"
        fields[2] = fields[2] or "."
        fields[3] = fields[3] or "N"
        fields[4] = fields[4] or "<X>"
        fields[5] = fields[5] or "."
        fields[6] = fields[6] or "."
        fields[7] = fields[7] or "."
        fields[8] = fields[8] or "GT"

        for idx in range(9, expected_cols):
            fields[idx] = fields[idx] or "./."

        fout.write("\t".join(fields) + "\n")

with open(meta_path, "w", encoding="utf-8") as meta:
    meta.write(f"{len(header_fields or [])}\n")
    meta.write(f"{total_rows}\n")
    meta.write(f"{patched_rows}\n")
PY

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
