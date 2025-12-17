#!/bin/bash
# =============================================================================
# Step 1B consolidated VCF validator
# Ensures each ChrXX_consolidated.vcf.gz is scheduler-ready for Step 1C
# =============================================================================
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: validate_consolidated_vcf.sh <vcf1.vcf.gz> [vcf2.vcf.gz ...]

Runs bcftools validation to catch malformed records (short rows, whitespace
issues, empty FORMAT/genotypes) immediately after GenotypeGVCFs produces
consolidated chromosome VCFs. Step 1C now performs its own normalization, so
this check simply ensures each file is readable end-to-end by bcftools.
EOF
}

if [ $# -lt 1 ]; then
    usage
    exit 1
fi

PIPELINE_ROOT_RESOLVED="${PIPELINE_ROOT:-}"
if [ -z "${PIPELINE_ROOT_RESOLVED}" ]; then
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    PIPELINE_ROOT_RESOLVED="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
fi
export PIPELINE_ROOT="${PIPELINE_ROOT_RESOLVED}"

BCFTOOLS_BIN="${BCFTOOLS_BIN:-bcftools}"

if command -v module >/dev/null 2>&1; then
    module purge
fi

if ! command -v "${BCFTOOLS_BIN}" >/dev/null 2>&1; then
    if command -v module >/dev/null 2>&1; then
        module load bcftools/1.18-GCC-12.3.0 >/dev/null 2>&1 || true
    fi
fi

if ! command -v "${BCFTOOLS_BIN}" >/dev/null 2>&1; then
    echo "[FATAL] bcftools not found; cannot validate consolidated VCFs." >&2
    exit 1
fi

status=0

for vcf in "$@"; do
    if [ ! -f "${vcf}" ]; then
        echo "[ERROR] VCF not found: ${vcf}" >&2
        status=1
        continue
    fi

    echo "[INFO] Validating ${vcf} with bcftools"
    if ! "${BCFTOOLS_BIN}" view -Ov -o /dev/null "${vcf}" >/dev/null 2>&1; then
        echo "[ERROR] bcftools reported an issue with ${vcf}" >&2
        status=1
        continue
    fi

done

exit "${status}"
