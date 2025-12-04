# Solution Plan for Persistent Beagle ninthTabPos Failure

## Problem
- Validator passes all 17 chromosomes
- Beagle crashes immediately on Chr01 with `ninthTabPos` error during header parsing (`VcfHeader.isDiploid()`)
- This indicates a **header-level** issue, not a data-row issue

## Root Cause Hypothesis
Beagle's `ninthTabPos()` literally counts tab characters to find the 9th tab (which separates INFO from FORMAT in a VCF). If the header or any part of the file has inconsistent whitespace, Beagle will fail even if field-based parsing works.

## Investigation Steps (Run on HPC)

### 1. Check Persisted Artifacts
```bash
cd /scratch/project/bigdata_apple/logs/test/step1c_debug/20251204_140918/
ls -lah
cat *.validate.log  # Should be empty or show "PASSED"
```

### 2. Inspect the Actual Chr01 File That Failed
```bash
# The file that Beagle rejected
zcat /scratch/temp/18482191/Chr01_consolidated.fixed.filtered.renamed.hdr.final.clean.sorted.beagle.vcf.gz | head -100 > ~/chr01_failed.txt

# Check for tab count in header
zcat /scratch/temp/18482191/Chr01_consolidated.fixed.filtered.renamed.hdr.final.clean.sorted.beagle.vcf.gz | grep "^#CHROM" | od -c | head -20

# Count tabs in #CHROM line
zcat /scratch/temp/18482191/Chr01_consolidated.fixed.filtered.renamed.hdr.final.clean.sorted.beagle.vcf.gz | grep "^#CHROM" | tr -cd '\t' | wc -c

# Check first 5 data lines
zcat /scratch/temp/18482191/Chr01_consolidated.fixed.filtered.renamed.hdr.final.clean.sorted.beagle.vcf.gz | grep -v "^#" | head -5 | od -c
```

### 3. Compare Working vs. Non-Working Chromosome
If other chromosomes work, compare:
```bash
# If Chr02 works but Chr01 doesn't
zcat /scratch/temp/18482191/Chr02*.beagle.vcf.gz | grep "^#CHROM" > ~/chr02_header.txt
zcat /scratch/temp/18482191/Chr01*.beagle.vcf.gz | grep "^#CHROM" > ~/chr01_header.txt
diff ~/chr01_header.txt ~/chr02_header.txt
```

## Quick Fix Options

### Option A: Strengthen the Validator (Most Likely)
Our validator checks `len(fields)` after split, but Beagle counts raw tab characters. We need to add:

1. **Strict tab count validation in header**
2. **Check that every data line has exactly expected_tabs tab characters**
3. **Reject files with ANY tab/space ambiguity**

### Option B: Bypass bcftools sort (If That's Corrupting)
The last step before Beagle is `bcftools sort`. If this is corrupting the VCF:
```bash
# Test without sorting
bcftools view -h your.vcf.gz  # Just check if basic structure is ok
```

### Option C: Use Different Compression
Some VCF tools create bgzip files that other tools can't read properly. Try:
```bash
# Decompress and recompress with standard gzip compression level
zcat file.vcf.gz | bgzip -c > file.recompressed.vcf.gz
```

## Recommended Implementation

### 1. Enhance the Validator (Immediate)
Add a strict byte-level check to `debug_beagle_vcf.sh`:

```python
# After reading header_line, check raw tab count
raw_header = header_line  # Before any processing
actual_tabs = raw_header.count('\t')
expected_tabs = len(header_fields) - 1

if actual_tabs != expected_tabs:
    print(f"[FAILED] Header has {actual_tabs} tabs but {expected_tabs} expected", file=sys.stderr)
    print(f"  This will cause Beagle ninthTabPos error", file=sys.stderr)
    sys.exit(1)

# For each data line, check tab count BEFORE split
for raw in fh:
    line = raw.rstrip("\n")
    tab_count = line.count('\t')
    if tab_count != expected_tabs:
        log_error(f"Line {line_count}: has {tab_count} tabs (expected {expected_tabs})")
```

### 2. Debug the Actual File on HPC
Since TMPDIR files are gone, we need to either:
- Re-run with `STEP1C_DEBUG_BEAGLE=true` to capture the failing files
- Or manually recreate the issue with a smaller test dataset

### 3. Add Pre-Beagle Sanity Check
Before launching Beagle, run a simple Java-based check or use Beagle's validation mode (if it has one).

## Testing Protocol

1. **Implement enhanced validator** (tab-count strict mode)
2. **Re-run step1c_submit.sh on test dataset**
3. **Check if new validator catches the issue**
4. **If validator still passes, manually inspect the Chr01 file** using investigation steps above
5. **Once we identify the exact malformation, add it to test fixtures**

## Next Actions

### For You (User)
Run the investigation steps above on the HPC to capture:
1. The actual `#CHROM` header line from the failing Chr01 file
2. The first few data lines
3. The exact tab counts

Provide this output and I'll pinpoint the exact issue.

### For Me (AI)
1. Update `debug_beagle_vcf.sh` to check raw tab counts (not just field counts after split)
2. Add this check to the validator immediately
3. Update documentation with this specific failure mode

