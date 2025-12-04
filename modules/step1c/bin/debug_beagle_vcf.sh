#!/bin/bash
# =============================================================================
# Beagle-focused VCF validator
# Detects issues that commonly trigger "VCF record format error (ninthTabPos)"
# =============================================================================
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: debug_beagle_vcf.sh <input.vcf.gz> [--no-snippets]

Validates a VCF for Beagle compatibility by checking:
  1. Presence and integrity of the #CHROM header
  2. Strict tab-delimited field counts per data row
  3. Non-empty CHROM/POS/REF/ALT/FORMAT columns
  4. FORMAT token syntax and genotype column counts
  5. Hidden whitespace (spaces or carriage returns) within records

Exit codes:
  0 - VCF passed all checks
  1 - VCF failed at least one check
  2 - Usage error (missing file, bad arguments)
EOF
}

if [ $# -lt 1 ]; then
    usage
    exit 2
fi

INPUT="$1"
INCLUDE_SNIPPETS=true
if [ "${2:-}" = "--no-snippets" ]; then
    INCLUDE_SNIPPETS=false
fi

if [ ! -f "${INPUT}" ]; then
    echo "[FATAL] Input VCF not found: ${INPUT}" >&2
    exit 2
fi

REPORT_FILE="$(mktemp)"
trap 'rm -f "${REPORT_FILE}"' EXIT

set +e
python3 - "$INPUT" "$INCLUDE_SNIPPETS" "$REPORT_FILE" <<'PY'
import gzip
import re
import sys

input_path = sys.argv[1]
include_snippets = sys.argv[2] == "True"
report_path = sys.argv[3]

errors = []
warnings = []

def log_error(msg, line=None):
    if line is not None and include_snippets:
        errors.append(f"{msg}\n  Line snippet: {line.rstrip()[:200]}")
    else:
        errors.append(msg)

def log_warn(msg, line=None):
    if line is not None and include_snippets:
        warnings.append(f"{msg}\n  Line snippet: {line.rstrip()[:200]}")
    else:
        warnings.append(msg)

try:
    fh = gzip.open(input_path, "rt", encoding="utf-8", errors="replace")
except OSError as exc:
    print(f"[FATAL] Unable to read {input_path}: {exc}", file=sys.stderr)
    sys.exit(2)

header_line = None
for raw in fh:
    if raw.startswith("#CHROM"):
        header_line = raw.rstrip("\n")
        break

if header_line is None:
    print("[FATAL] #CHROM header not found", file=sys.stderr)
    sys.exit(1)

header_fields = header_line.split("\t")
tab_count = header_line.count("\t")

if len(header_fields) < 9:
    print(f"[FAILED] Header has only {len(header_fields)} columns (need >=9)", file=sys.stderr)
    sys.exit(1)

sample_count = len(header_fields) - 9
expected_cols = 9 + sample_count
expected_tabs = expected_cols - 1

print(f"Input: {input_path}")
print(f"Header columns: {len(header_fields)} (tabs: {tab_count}, expected tabs: {expected_tabs})")
print(f"Sample count: {sample_count}")

if tab_count != expected_tabs:
    log_error(f"Header has {tab_count} tab characters but {expected_tabs} expected")
    log_error(f"Beagle will fail with ninthTabPos error on this file")
    print(f"[FAILED] Header tab count mismatch: {tab_count} != {expected_tabs}", file=sys.stderr)
    sys.exit(1)

line_count = 0

for raw in fh:
    if raw.startswith("#"):
        continue

    line = raw.rstrip("\n")
    line_count += 1

    if "\r" in line:
        log_warn(f"Line {line_count}: contains carriage return (\\r)", line)

    if re.search(r"[^\t] [^\t]", line):
        log_warn(f"Line {line_count}: contains spaces between non-tab characters (mixed delimiters)", line)

    # CRITICAL: Check raw tab count before parsing (this is what Beagle's ninthTabPos checks)
    actual_tabs = line.count('\t')
    if actual_tabs != expected_tabs:
        log_error(f"Line {line_count}: has {actual_tabs} tabs (expected {expected_tabs}) - Beagle ninthTabPos will fail", line)

    fields = line.split("\t")

    if len(fields) != expected_cols:
        log_error(f"Line {line_count}: has {len(fields)} columns (expected {expected_cols})", line)
        continue

    if not fields[0] or not fields[1] or not fields[3] or not fields[4]:
        log_error(f"Line {line_count}: empty CHROM/POS/REF/ALT field", line)

    format_field = fields[8]
    if not format_field:
        log_error(f"Line {line_count}: empty FORMAT column", line)
    elif not re.fullmatch(r"[A-Za-z0-9_:.+-]+", format_field):
        log_error(f"Line {line_count}: invalid FORMAT token '{format_field}'", line)

    for idx, sample in enumerate(fields[9:], start=9):
        if sample == "":
            log_error(f"Line {line_count}: sample column {idx+1} is empty", line)
            break
        if sample == ".":
            continue
        if not re.match(r"[0-9./|]+(:.*)?$", sample):
            log_warn(f"Line {line_count}: unusual genotype token '{sample}' in column {idx+1}", line)

fh.close()

with open(report_path, "w", encoding="utf-8") as rep:
    if errors:
        rep.write("Errors:\n-------\n")
        rep.write("\n".join(errors[:20]))
        if len(errors) > 20:
            rep.write(f"\n... {len(errors) - 20} more errors not shown\n")
        rep.write("\n\n")
    if warnings:
        rep.write("Warnings:\n---------\n")
        rep.write("\n".join(warnings[:20]))
        if len(warnings) > 20:
            rep.write(f"\n... {len(warnings) - 20} more warnings not shown\n")
        rep.write("\n")

print(f"Validated data lines: {line_count}")

if errors:
    print("[FAILED] VCF is not Beagle-compatible")
    sys.exit(1)

if warnings:
    print("[WARN] VCF contains potential issues; see report for details")

print("[PASSED] VCF passed Beagle validation")
PY
validator_status=$?
set -e

if [ -s "${REPORT_FILE}" ]; then
    echo "------------------------------------------"
    cat "${REPORT_FILE}"
    echo "------------------------------------------"
fi

exit "${validator_status}"

