# Immediate Diagnostic Steps for ninthTabPos Failure

## The Situation
- Validator passes all chromosomes ✓
- Tab counts are correct ✓
- Beagle still fails on Chr01 immediately ✗

**This means:** The issue is NOT tab count - it's something else Beagle checks in `VcfHeader.isDiploid()`.

## Critical: Check Persisted Validator Logs First

```bash
# Go to the step1c_debug directory from the latest run
cd /scratch/project/bigdata_apple/logs/test/step1c_debug/
ls -lt | head -5

# Look at the most recent run (should be 20251204_142xxx based on timestamps)
cd 20251204_142240  # Or whatever the actual timestamp is

# Check what was captured
ls -lah

# Read the validator log for Chr01
cat Chr01_consolidated.fixed.filtered.renamed.hdr.final.clean.sorted.beagle.validate.log
```

**Expected:** The validator log should show "PASSED" but we need to see the actual output.

## If Validator Log Says "PASSED": Inspect the Actual File

The issue is likely in how `bcftools sort` is writing the file. Let's examine the raw bytes:

### 1. Check the Header Line Structure
```bash
# Look at JUST the #CHROM line (raw bytes)
zcat /scratch/temp/18482343/Chr01_consolidated.fixed.filtered.renamed.hdr.final.clean.sorted.beagle.vcf.gz | \
  grep "^#CHROM" | od -An -tx1 | head -5

# This shows hex bytes - look for:
# 09 = tab character (should appear exactly 10 times for 11 columns)
# 20 = space character (should NOT appear between sample names)
```

### 2. Compare Before and After bcftools sort
```bash
# Check the file BEFORE bcftools sort
zcat /scratch/temp/18482343/Chr01_consolidated.fixed.filtered.renamed.hdr.final.clean.vcf.gz | \
  grep "^#CHROM" | od -An -tx1 | head -5

# Check the file AFTER bcftools sort
zcat /scratch/temp/18482343/Chr01_consolidated.fixed.filtered.renamed.hdr.final.clean.sorted.beagle.vcf.gz | \
  grep "^#CHROM" | od -An -tx1 | head -5

# If these differ, bcftools sort is corrupting the header
```

### 3. Test If It's a Compression Issue
```bash
# Try decompressing and recompressing
cd /scratch/temp/18482343/
zcat Chr01_consolidated.fixed.filtered.renamed.hdr.final.clean.sorted.beagle.vcf.gz | \
  bgzip -c > Chr01_test.vcf.gz
tabix -f Chr01_test.vcf.gz

# Try running Beagle on the recompressed file
java -jar $BEAGLE_JAR gt=Chr01_test.vcf.gz out=chr01_test_output
```

### 4. Look at Sample Names in Header
```bash
# Extract just the sample names
zcat /scratch/temp/18482343/Chr01_consolidated.fixed.filtered.renamed.hdr.final.clean.sorted.beagle.vcf.gz | \
  grep "^#CHROM" | cut -f10-

# Check for any special characters, spaces, or unusual encodings
zcat /scratch/temp/18482343/Chr01_consolidated.fixed.filtered.renamed.hdr.final.clean.sorted.beagle.vcf.gz | \
  grep "^#CHROM" | cut -f10- | od -An -tx1
```

### 5. Check First Data Line
```bash
# The error happens in isDiploid() which checks the first GT field
zcat /scratch/temp/18482343/Chr01_consolidated.fixed.filtered.renamed.hdr.final.clean.sorted.beagle.vcf.gz | \
  grep -v "^#" | head -1 | od -An -tx1 | head -10
```

## Known Beagle ninthTabPos Triggers

Based on Beagle source code, `ninthTabPos()` fails when:

1. **Line has fewer than 9 tabs** (we already check this)
2. **Header line structure doesn't match data line structure**
3. **Sample names contain spaces or tabs**
4. **Malformed UTF-8 or special characters**
5. **Inconsistent line endings** (mixed \n and \r\n)

## Quick Test: Try a Working Chromosome

If Chr02 or another chromosome works:

```bash
# Compare Chr01 vs Chr02 headers byte-by-byte
zcat /scratch/temp/18482343/Chr01_*.beagle.vcf.gz | grep "^#CHROM" > /tmp/chr01_header.txt
zcat /scratch/temp/18482343/Chr02_*.beagle.vcf.gz | grep "^#CHROM" > /tmp/chr02_header.txt

# Should be identical (just different sample IDs if any)
diff /tmp/chr01_header.txt /tmp/chr02_header.txt

# If different, show me the difference
diff -u /tmp/chr01_header.txt /tmp/chr02_header.txt | od -An -tx1
```

## Nuclear Option: Skip bcftools sort

The last transformation before Beagle is `bcftools sort`. Let's test without it:

```bash
# In step1c_job.sh, comment out the bcftools sort step and run:
# BEAGLE_VCF_PATHS+=("${sorted_path%.beagle.vcf.gz}.vcf.gz")
# Instead of:
# BEAGLE_VCF_PATHS+=("${sorted_path}")
```

Or manually test now:

```bash
cd /scratch/temp/18482343/
# Use the file BEFORE bcftools sort
java -Xmx8g -jar $BEAGLE_JAR \
  gt=Chr01_consolidated.fixed.filtered.renamed.hdr.final.clean.vcf.gz \
  out=chr01_unsorted_test
```

## Report Back To Me

Please run these commands and send me:
1. The output of the validator log (should be in step1c_debug directory)
2. The hex dump of the #CHROM line
3. Whether Chr02 or other chromosomes have the same issue
4. Whether the unsorted VCF works with Beagle

This will tell us exactly what's wrong.

