#!/bin/bash
# =============================================================================
# Step 1C VCF repair + validator smoke test
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
FIX_SCRIPT="${PIPELINE_DIR}/modules/step1c/bin/fix_vcf_fill_missing.sh"
VALIDATOR="${PIPELINE_DIR}/modules/step1c/bin/debug_beagle_vcf.sh"
FIXTURES=(
    "${PIPELINE_DIR}/test/fixtures/malformed_chr01.vcf"
    "${PIPELINE_DIR}/test/fixtures/malformed_chr01_empty_sample.vcf"
)

for tool in bgzip tabix python3; do
    if ! command -v "${tool}" >/dev/null 2>&1; then
        echo "[SKIP] ${tool} not available; skipping VCF repair test."
        exit 0
    fi
done

if [ ! -x "${FIX_SCRIPT}" ] || [ ! -x "${VALIDATOR}" ]; then
    echo "[ERROR] Required scripts missing or not executable."
    exit 1
fi

TMP_DIR="$(mktemp -d)"
trap 'rm -rf "${TMP_DIR}"' EXIT

for fixture in "${FIXTURES[@]}"; do
    if [ ! -f "${fixture}" ]; then
        echo "[ERROR] Fixture not found: ${fixture}"
        exit 1
    fi

    base_name="$(basename "${fixture}" .vcf)"
    raw_gz="${TMP_DIR}/${base_name}.vcf.gz"

    python3 - <<'PY'
import gzip, sys, pathlib
src = pathlib.Path(sys.argv[1])
dst = pathlib.Path(sys.argv[2])
with gzip.open(dst, "wt") as out, open(src, "r", encoding="utf-8") as fh:
    out.write(fh.read())
PY "${fixture}" "${raw_gz}"

    echo "[TEST] Expect validator to fail on ${base_name}..."
    if "${VALIDATOR}" "${raw_gz}" > "${TMP_DIR}/${base_name}.raw.log" 2>&1; then
        echo "[ERROR] Validator unexpectedly passed malformed input (${base_name})."
        cat "${TMP_DIR}/${base_name}.raw.log"
        exit 1
    fi

    fixed_gz="${TMP_DIR}/${base_name}.fixed.vcf.gz"
    echo "[TEST] Running fix_vcf_fill_missing.sh on ${base_name}..."
    if ! "${FIX_SCRIPT}" "${raw_gz}" "${fixed_gz}" > "${TMP_DIR}/${base_name}.fix.log" 2>&1; then
        echo "[ERROR] fix_vcf_fill_missing.sh failed (${base_name}):"
        cat "${TMP_DIR}/${base_name}.fix.log"
        exit 1
    fi

    echo "[TEST] Validator should now pass for ${base_name}..."
    if ! "${VALIDATOR}" "${fixed_gz}" > "${TMP_DIR}/${base_name}.fixed.log" 2>&1; then
        echo "[ERROR] Validator failed on repaired VCF (${base_name}):"
        cat "${TMP_DIR}/${base_name}.fixed.log"
        exit 1
    fi
done

echo "[PASS] VCF repair + validator pipeline succeeded for all fixtures."

