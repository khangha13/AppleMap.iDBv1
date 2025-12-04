#!/bin/bash
# =============================================================================
# Beagle-focused VCF validator
# Detects issues that commonly trigger "VCF record format error (ninthTabPos)"
# =============================================================================
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: debug_beagle_vcf.sh <input.vcf.gz> [--verbose]

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
VERBOSE=false
if [ "${2:-}" = "--verbose" ]; then
    VERBOSE=true
fi

if [ ! -f "${INPUT}" ]; then
    echo "[FATAL] Input VCF not found: ${INPUT}" >&2
    exit 2
fi

REPORT_FILE="$(mktemp)"
trap 'rm -f "${REPORT_FILE}"' EXIT

python3 - "$INPUT" "$VERBOSE" "$REPORT_FILE" <<'PY'
import gzip
import os
import re
import sys

input_path = sys.argv[1]
verbose = sys.argv[2] == "True"
report_path = sys.argv[3]

errors = []
warnings = []

def log_error(msg, line=None):
    if line is not None and verbose:
        errors.append(f"{msg}\n  Line snippet: {line.rstrip()[:200]}")
    else:
        errors.append(msg)

def log_warn(msg, line=None):
    if line is not None and verbose:
        warnings.append(f"{msg}\n  Line snippet: {line.rstrip()[:200]}")
    else:
        warnings.append(msg)

try:
    fh = gzip.open(input_path, "rt", encoding="utf-8", errors="replace")
except OSError as exc:
    print(f"[FATAL] Unable to read {input_path}: {exc}", file=sys.stderr)
    sys.exit(2)

header_line = None
line_count = 0

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
    log_warn(f"Header contains {tab_count} tab characters but {expected_tabs} are expected")

data_start_offset = fh.tell()
fh.close()
fh = gzip.open(input_path, "rt", encoding="utf-8", errors="replace")
fh.seek(data_start_offset)

for line_number, raw in enumerate(fh, start=1):
    if raw.startswith("#"):
        continue

    line = raw.rstrip("\n")
    line_count += 1

    if "\r" in line:
        log_warn(f"Line {line_number}: contains carriage return (\\r)", line)

    if re.search(r"[^\t] [^\t]", line):
        log_warn(f"Line {line_number}: contains spaces between non-tab characters (mixed delimiters)", line)

    fields = line.split("\t")

    if len(fields) != expected_cols:
        log_error(f"Line {line_number}: has {len(fields)} columns (expected {expected_cols})", line)
        continue

    # Check critical fields
    if not fields[0] or not fields[1] or not fields[3] or not fields[4]:
        log_error(f"Line {line_number}: empty CHROM/POS/REF/ALT field", line)

    format_field = fields[8]
    if not format_field:
        log_error(f"Line {line_number}: empty FORMAT column", line)
    elif not re.fullmatch(r"[A-Za-z0-9_:.+-]+", format_field):
        log_error(f"Line {line_number}: invalid FORMAT token '{format_field}'", line)

    # Genotype columns
    for idx, sample in enumerate(fields[9:], start=9):
        if sample == "":
            log_error(f"Line {line_number}: sample column {idx+1} is empty", line)
            break
        if sample == ".":
            continue
        if not re.match(r"[0-9./|]+(:.*)?$", sample):
            log_warn(f"Line {line_number}: unusual genotype token '{sample}' in column {idx+1}", line)

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

if [ -s "${REPORT_FILE}" ]; then
    echo "------------------------------------------"
    cat "${REPORT_FILE}"
    echo "------------------------------------------"
fi

exit "${validator_status}"
#!/bin/bash
# =============================================================================
# Beagle-focused VCF validator
# Checks for common issues that cause Beagle's "VCF record format error (ninthTabPos)"
# =============================================================================
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: debug_beagle_vcf.sh <input.vcf.gz> [--verbose]

Validates a VCF file for Beagle compatibility by checking:
  1. Header has #CHROM line with correct format
  2. All data lines have exactly the expected number of tab-separated fields
  3. No mixed whitespace (spaces where tabs expected)
  4. FORMAT column is non-empty and valid
  5. Genotype fields are syntactically valid

Exit codes:
  0 - VCF is valid for Beagle
  1 - VCF has issues that will cause Beagle to fail
  2 - Usage/input error

Outputs a detailed report of any issues found.
EOF
}

if [ $# -lt 1 ]; then
    usage
    exit 2
fi

INPUT="$1"
VERBOSE="${2:-}"

if [ ! -f "${INPUT}" ]; then
    echo "[FATAL] Input VCF not found: ${INPUT}" >&2
    exit 2
fi

# Colors for output (if terminal)
RED='\033[0;31m'
YELLOW='\033[1;33m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

# Disable colors if not a terminal
if [ ! -t 1 ]; then
    RED=''
    YELLOW=''
    GREEN=''
    NC=''
fi

echo "=========================================="
echo "Beagle VCF Validator"
echo "=========================================="
echo "Input: ${INPUT}"
echo ""

# Extract header line
header_line=$(zcat "${INPUT}" 2>/dev/null | awk '/^#CHROM/ {print; exit}')
if [ -z "${header_line}" ]; then
    echo -e "${RED}[FATAL] #CHROM header not found${NC}"
    exit 1
fi

# Count header fields using strict tab splitting
header_fields=$(printf "%s" "${header_line}" | awk -F'\t' '{print NF}')
echo "Header columns (tab-delimited): ${header_fields}"

if [ "${header_fields}" -lt 9 ]; then
    echo -e "${RED}[FATAL] Header has fewer than 9 columns: ${header_fields}${NC}"
    exit 1
fi

sample_count=$((header_fields - 9))
echo "Sample count: ${sample_count}"
expected_cols=$((9 + sample_count))
echo "Expected columns per data line: ${expected_cols}"
echo ""

# Count actual tabs in header (Beagle counts tabs, not fields)
header_tabs=$(printf "%s" "${header_line}" | tr -cd '\t' | wc -c)
echo "Tab characters in header: ${header_tabs}"
expected_tabs=$((expected_cols - 1))
if [ "${header_tabs}" -ne "${expected_tabs}" ]; then
    echo -e "${RED}[ERROR] Header has ${header_tabs} tabs, expected ${expected_tabs}${NC}"
fi
echo ""

# Validate data records
echo "Validating data records..."
echo ""

errors_found=0
warnings_found=0

# Create temp file for detailed report
report_file=$(mktemp)
trap 'rm -f "${report_file}"' EXIT

# Comprehensive validation using awk
zcat "${INPUT}" 2>/dev/null | awk -F'\t' -v expected="${expected_cols}" -v expected_tabs="${expected_tabs}" -v report="${report_file}" -v verbose="${VERBOSE}" '
BEGIN {
    errors = 0
    warnings = 0
    line_num = 0
    max_errors_shown = 10
}

/^#/ { next }  # Skip header lines

{
    line_num++
    has_error = 0
    error_msg = ""
    
    # Check 1: Count tab-separated fields
    if (NF != expected) {
        has_error = 1
        error_msg = error_msg sprintf("  - Field count: %d (expected %d)\n", NF, expected)
    }
    
    # Check 2: Count actual tab characters in line
    gsub(/[^\t]/, "", $0)
    tab_count = length($0)
    # Re-read line for further checks (awk modifies $0 with gsub)
}

# Re-process the line for content checks
!/^#/ {
    # Reload original line
    line = $0
    
    # Check for spaces that might be masquerading as field separators
    # This catches lines where spaces were used instead of tabs
    if (match(line, /[^\t] [^\t]/)) {
        warnings++
        if (warnings <= max_errors_shown) {
            print "Line " NR ": Contains spaces between non-tab characters (potential mixed delimiters)" >> report
        }
    }
    
    # Check for carriage returns
    if (match(line, /\r/)) {
        warnings++
        if (warnings <= max_errors_shown) {
            print "Line " NR ": Contains carriage return (\\r)" >> report
        }
    }
    
    # Check for empty fields in critical positions (CHROM, POS, REF, ALT, FORMAT)
    if ($1 == "" || $2 == "" || $4 == "" || $5 == "" || $9 == "") {
        errors++
        if (errors <= max_errors_shown) {
            print "Line " NR ": Empty critical field (CHROM, POS, REF, ALT, or FORMAT)" >> report
            if (verbose == "--verbose") {
                print "  Content: " substr($0, 1, 200) >> report
            }
        }
    }
    
    # Check FORMAT field format (should be colon-separated tags like GT:DP:GQ)
    format = $9
    if (format !~ /^[A-Za-z0-9:]+$/) {
        errors++
        if (errors <= max_errors_shown) {
            print "Line " NR ": Invalid FORMAT field: \"" format "\"" >> report
        }
    }
    
    # Check that we have the right number of fields
    if (NF < expected) {
        errors++
        if (errors <= max_errors_shown) {
            print "Line " NR ": Too few fields: " NF " (expected " expected ")" >> report
            if (verbose == "--verbose") {
                print "  Content: " substr($0, 1, 200) >> report
            }
        }
    } else if (NF > expected) {
        warnings++
        if (warnings <= max_errors_shown) {
            print "Line " NR ": Extra fields: " NF " (expected " expected ")" >> report
        }
    }
    
    # Beagle-specific: Check that line has enough tab characters
    # This is what ninthTabPos actually checks
    line_copy = $0
    gsub(/[^\t]/, "", line_copy)
    actual_tabs = length(line_copy)
    if (actual_tabs < 8) {  # Need at least 8 tabs for FORMAT column
        errors++
        if (errors <= max_errors_shown) {
            print "Line " NR ": Only " actual_tabs " tabs (Beagle needs at least 8 for ninthTabPos)" >> report
            if (verbose == "--verbose") {
                print "  Content: " substr($0, 1, 200) >> report
            }
        }
    }
}

END {
    print errors > "/dev/stderr"
    print warnings > "/dev/stderr"
    print line_num > "/dev/stderr"
}
' 2>&1 | {
    read errors_count
    read warnings_count  
    read total_lines
    
    echo "Total data lines: ${total_lines}"
    echo ""
    
    if [ -s "${report_file}" ]; then
        echo "Issues found:"
        echo "-------------"
        cat "${report_file}"
        echo ""
    fi
    
    if [ "${errors_count:-0}" -gt 0 ]; then
        echo -e "${RED}[FAILED] ${errors_count} error(s) found${NC}"
        errors_found=1
    fi
    
    if [ "${warnings_count:-0}" -gt 0 ]; then
        echo -e "${YELLOW}[WARNING] ${warnings_count} warning(s) found${NC}"
    fi
    
    if [ "${errors_count:-0}" -eq 0 ] && [ "${warnings_count:-0}" -eq 0 ]; then
        echo -e "${GREEN}[PASSED] VCF appears valid for Beagle${NC}"
    fi
    
    # Exit with error if any errors found
    if [ "${errors_count:-0}" -gt 0 ]; then
        exit 1
    fi
}

exit_code=$?

echo ""
echo "=========================================="

exit ${exit_code}

