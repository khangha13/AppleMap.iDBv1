# Detailed Explanation: Filtering Options for Biallelic VCFs

## Overview

When passing a VCF file that is already biallelic SNPs into the `master_vcf_analysis.sh` pipeline, the current implementation will redundantly filter it again using `bcftools view -m2 -M2 -v snps`. This document explains two options to avoid this redundant filtering.

**Note (2025-12-19):** Step 1D now always filters to biallelic SNPs for consistency. The options below are retained as historical/reference material.

---

## Option 1: Check if Input VCF is Already Biallelic SNPs

### How It Works

This option performs a **runtime check** on the actual VCF file content to determine if it's already filtered to biallelic SNPs before running the filtering step.

### Implementation Details

#### Step 1: Sample Variants from the VCF
- Use `bcftools query` to extract a sample of variants (e.g., first 1000 lines)
- Check three properties:
  1. **Allele count**: Are all variants biallelic? (2 alleles total)
  2. **Variant type**: Are all variants SNPs? (not indels)
  3. **Completeness**: Do we have enough variants to make a reliable check?

#### Step 2: Check Variant Properties
```bash
# Example check command
bcftools query -f '%ALT\n' "${vcf_file}" | head -1000 | \
  awk '
    {
      split($0, alts, ",");
      if (length(alts) > 1) {  # More than 1 allele = multiallelic
        exit 1;
      }
    }'
```

#### Step 3: Check for Non-SNP Variants
```bash
# Check if any variants are indels (not SNPs)
bcftools query -f '%TYPE\n' "${vcf_file}" | head -1000 | \
  awk '/INDEL/ { exit 1 }'
```

#### Step 4: Decision Logic
```bash
is_already_biallelic_snps() {
    local vcf="$1"
    local sample_size=1000
    
    # Check 1: Sample variants and count alleles
    local multiallelic_count=$(bcftools query -f '%ALT\n' "${vcf}" 2>/dev/null | \
        head -${sample_size} | \
        awk '{split($0, alts, ","); if (length(alts) > 1) count++} END {print count+0}')
    
    # Check 2: Check for indels
    local indel_count=$(bcftools query -f '%TYPE\n' "${vcf}" 2>/dev/null | \
        head -${sample_size} | \
        awk '/INDEL/ {count++} END {print count+0}')
    
    # If no multiallelic variants and no indels in sample, assume it's already filtered
    if [ "${multiallelic_count}" -eq 0 ] && [ "${indel_count}" -eq 0 ]; then
        return 0  # Already biallelic SNPs
    else
        return 1  # Needs filtering
    fi
}
```

#### Step 5: Conditional Filtering
```bash
if is_already_biallelic_snps "${vcf_file}"; then
    log_info "${chr_name}: VCF appears to be already filtered to biallelic SNPs, skipping filtering step"
    # Use input VCF directly
    cp "${vcf_file}" "${filtered_output}"
    bcftools index -t "${filtered_output}"
else
    log_info "${chr_name}: filtering to biallelic SNPs"
    bcftools view -m2 -M2 -v snps -Oz -o "${filtered_temp}" "${vcf_file}"
    # ... rest of filtering logic
fi
```

### Advantages
- ✅ **Accurate**: Checks actual file content, not just filename
- ✅ **Flexible**: Works regardless of filename convention
- ✅ **Safe**: Can handle edge cases where filename doesn't match content

### Disadvantages
- ⚠️ **Performance**: Adds overhead (reading VCF file twice - once for check, once for filtering)
- ⚠️ **Complexity**: More code to maintain
- ⚠️ **False positives**: If file has some multiallelic variants mixed in, might incorrectly skip filtering
- ⚠️ **Sampling**: Only checks a sample, not the entire file (though unlikely to have issues)

### Performance Impact
- **Time**: ~5-30 seconds per chromosome (depending on VCF size) for the check
- **I/O**: Reads ~1000 lines from VCF file
- **Benefit**: Saves filtering time if VCF is already filtered (~2-10 minutes per chromosome)

---

## Option 2: Skip Filtering Based on Filename Pattern

### How It Works

This option uses **filename heuristics** to determine if a VCF is likely already filtered before running the filtering step.

### Implementation Details

#### Step 1: Check Filename Patterns
Common patterns that suggest a VCF is already filtered:
- `*_snps.vcf.gz` - Contains "snps" in filename
- `*_biallelic.vcf.gz` - Contains "biallelic" in filename
- `*_filtered.vcf.gz` - Contains "filtered" in filename
- `*_snp.vcf.gz` - Contains "snp" in filename

#### Step 2: Pattern Matching Logic
```bash
is_likely_already_filtered() {
    local filename="$1"
    local basename=$(basename "${filename}" .vcf.gz)
    local basename_lower=$(echo "${basename}" | tr '[:upper:]' '[:lower:]')
    
    # List of patterns that suggest filtering
    local patterns=("snps" "snp" "biallelic" "filtered" "snv")
    
    for pattern in "${patterns[@]}"; do
        if [[ "${basename_lower}" == *"${pattern}"* ]]; then
            return 0  # Likely already filtered
        fi
    done
    
    return 1  # Probably needs filtering
}
```

#### Step 3: Conditional Filtering
```bash
if is_likely_already_filtered "${vcf_basename}"; then
    log_info "${chr_name}: Filename suggests VCF is already filtered (contains 'snps'/'biallelic'/'filtered'), skipping filtering step"
    # Use input VCF directly
    cp "${vcf_file}" "${filtered_output}"
    bcftools index -t "${filtered_output}"
else
    log_info "${chr_name}: filtering to biallelic SNPs"
    bcftools view -m2 -M2 -v snps -Oz -o "${filtered_temp}" "${vcf_file}"
    # ... rest of filtering logic
fi
```

#### Step 4: Optional: Add Command-Line Flag
You could also add a flag to explicitly skip filtering:
```bash
--skip-filtering    Skip biallelic SNP filtering step (assumes input is already filtered)
```

### Advantages
- ✅ **Fast**: No file I/O, just string matching
- ✅ **Simple**: Easy to implement and understand
- ✅ **Predictable**: User controls behavior via filename
- ✅ **Zero overhead**: No performance penalty

### Disadvantages
- ⚠️ **Not guaranteed**: Filename doesn't guarantee file content
- ⚠️ **Brittle**: Requires consistent naming conventions
- ⚠️ **User-dependent**: Relies on users naming files correctly
- ⚠️ **False negatives**: If user doesn't follow naming convention, will filter unnecessarily

### Performance Impact
- **Time**: < 1 millisecond (negligible)
- **I/O**: None
- **Benefit**: Saves filtering time if filename matches pattern (~2-10 minutes per chromosome)

---

## Comparison Table

| Aspect | Option 1 (Content Check) | Option 2 (Filename Pattern) |
|--------|-------------------------|----------------------------|
| **Accuracy** | High (checks actual content) | Medium (relies on naming) |
| **Speed** | Slow (~5-30 sec per chrom) | Fast (< 1 ms) |
| **Complexity** | Higher | Lower |
| **Reliability** | High (works regardless of name) | Medium (depends on naming) |
| **False Positives** | Low (rare) | Medium (possible) |
| **False Negatives** | Low (rare) | Medium (if naming not followed) |
| **Maintenance** | More code to maintain | Less code |

---

## Recommended Approach: Hybrid Option

### Best of Both Worlds

Combine both approaches:
1. **First check filename** (fast, Option 2)
2. **If filename doesn't match, check content** (Option 1)
3. **Add explicit flag** (`--skip-filtering`) for user control

### Implementation Pseudocode
```bash
skip_filtering=false

# Check 1: Command-line flag
if [ "${SKIP_FILTERING}" = "true" ]; then
    skip_filtering=true
fi

# Check 2: Filename pattern (fast)
if [ "${skip_filtering}" = false ] && is_likely_already_filtered "${vcf_basename}"; then
    log_info "${chr_name}: Filename suggests already filtered, skipping filter step"
    skip_filtering=true
fi

# Check 3: Content check (only if filename didn't match)
if [ "${skip_filtering}" = false ]; then
    if is_already_biallelic_snps "${vcf_file}"; then
        log_info "${chr_name}: Content check confirms already filtered, skipping filter step"
        skip_filtering=true
    fi
fi

# Execute filtering or skip
if [ "${skip_filtering}" = true ]; then
    # Use input directly
    cp "${vcf_file}" "${filtered_output}"
    bcftools index -t "${filtered_output}"
else
    # Run filtering
    bcftools view -m2 -M2 -v snps -Oz -o "${filtered_temp}" "${vcf_file}"
    # ... rest of logic
fi
```

---

## Example Usage Scenarios

### Scenario 1: VCF named `Chr00_snps.vcf.gz` (already filtered)
- **Option 1**: Checks content → Finds it's already filtered → Skips filtering ✅
- **Option 2**: Sees "snps" in filename → Skips filtering ✅
- **Hybrid**: Sees "snps" in filename → Skips filtering immediately ✅

### Scenario 2: VCF named `Chr00.vcf.gz` but content is already filtered
- **Option 1**: Checks content → Finds it's already filtered → Skips filtering ✅
- **Option 2**: No pattern match → Filters anyway (redundant) ❌
- **Hybrid**: No filename match → Checks content → Skips filtering ✅

### Scenario 3: VCF named `Chr00_snps.vcf.gz` but contains multiallelic variants
- **Option 1**: Checks content → Finds multiallelic → Filters anyway ✅
- **Option 2**: Sees "snps" in filename → Skips filtering (incorrect) ❌
- **Hybrid**: Sees "snps" → Skips (but could add content check as safety) ⚠️

---

## Recommendation

**Implement Option 2 (Filename Pattern) first** because:
1. It's simpler and faster
2. Most users will name filtered files appropriately
3. Low risk if implemented carefully
4. Can add Option 1 later if needed

**Consider adding `--skip-filtering` flag** for explicit user control when they know their VCFs are already filtered.
