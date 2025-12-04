#!/bin/bash
# =============================================================================
# Step 1C VCF repair + validator smoke test
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
FIX_SCRIPT="${PIPELINE_DIR}/modules/step1c/bin/fix_vcf_fill_missing.sh"
VALIDATOR="${PIPELINE_DIR}/modules/step1c/bin/debug_beagle_vcf.sh"
FIXTURE="${PIPELINE_DIR}/test/fixtures/malformed_chr01.vcf"

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

if [ ! -f "${FIXTURE}" ]; then
    echo "[ERROR] Fixture not found: ${FIXTURE}"
    exit 1
fi

TMP_DIR="$(mktemp -d)"
trap 'rm -rf "${TMP_DIR}"' EXIT

RAW_GZ="${TMP_DIR}/malformed.vcf.gz"
python3 - <<'PY'
import gzip, sys, pathlib
fixture = pathlib.Path(sys.argv[1])
output = pathlib.Path(sys.argv[2])
with gzip.open(output, "wt") as out, open(fixture, "r", encoding="utf-8") as src:
    out.write(src.read())
PY "${FIXTURE}" "${RAW_GZ}"

echo "[TEST] Expecting validator to fail on malformed VCF..."
if "${VALIDATOR}" "${RAW_GZ}" > "${TMP_DIR}/raw_validate.log" 2>&1; then
    echo "[ERROR] Validator unexpectedly passed malformed input."
    cat "${TMP_DIR}/raw_validate.log"
    exit 1
fi

FIXED_GZ="${TMP_DIR}/malformed.fixed.vcf.gz"
echo "[TEST] Running fix_vcf_fill_missing.sh..."
if ! "${FIX_SCRIPT}" "${RAW_GZ}" "${FIXED_GZ}" > "${TMP_DIR}/fix.log" 2>&1; then
    echo "[ERROR] fix_vcf_fill_missing.sh failed:"
    cat "${TMP_DIR}/fix.log"
    exit 1
fi

echo "[TEST] Validator should now pass..."
if ! "${VALIDATOR}" "${FIXED_GZ}" > "${TMP_DIR}/fixed_validate.log" 2>&1; then
    echo "[ERROR] Validator failed on repaired VCF:"
    cat "${TMP_DIR}/fixed_validate.log"
    exit 1
fi

echo "[PASS] VCF repair + validator pipeline succeeded."

