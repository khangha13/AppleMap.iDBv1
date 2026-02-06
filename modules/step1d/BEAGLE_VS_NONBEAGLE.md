# BEAGLE vs Non-BEAGLE Mode Comparison

**Purpose:** Comprehensive guide on the differences between BEAGLE and non-BEAGLE modes in Step1D

---

## Quick Summary

| Aspect | Non-BEAGLE (Default) | BEAGLE Mode |
|--------|----------------------|-------------|
| **Data Type** | Raw VCFs from GATK variant calling | Imputed VCFs from BEAGLE software |
| **Input Files** | 18 files (Chr00-Chr17) | 17 files (Chr01-Chr17, no Chr00) |
| **Key Metrics** | Depth (DP), Quality by Depth (QD), Mapping Quality (MQ) | Allele Frequency (AF), Dosage R² (DR2), Imputation flag (IMP) |
| **Plots Generated** | 14+ plot types (including depth-dependent) | 8 plot types (no depth-dependent plots) |
| **Use Case** | QC for sequencing data quality | QC for imputation quality |

---

## What is BEAGLE?

**BEAGLE** is a software tool for genotype imputation and phasing. It:
- Fills in missing genotypes using reference panels
- Improves low-coverage sequencing data
- Produces **imputed** VCFs with special INFO tags

**When to use BEAGLE mode in Step1D:**
- ✅ Your VCFs were processed through BEAGLE imputation
- ✅ You want to assess imputation quality (DR2 scores)
- ✅ Your VCFs have `AF`, `DR2`, `IMP` INFO tags

**When NOT to use BEAGLE mode:**
- ❌ Raw VCFs directly from GATK HaplotypeCaller
- ❌ VCFs from other variant callers without imputation
- ❌ You want depth-based QC metrics

---

## Detailed Comparison

### 1. Input Files Expected

#### Non-BEAGLE (Default)
```bash
Chr00.vcf.gz  # Chromosome 0 (scaffold/unplaced)
Chr01.vcf.gz  # Chromosome 1
Chr02.vcf.gz  # Chromosome 2
...
Chr17.vcf.gz  # Chromosome 17
```
**Total:** 18 VCF files

#### BEAGLE Mode
```bash
Chr01.vcf.gz  # Chromosome 1
Chr02.vcf.gz  # Chromosome 2
...
Chr17.vcf.gz  # Chromosome 17
```
**Total:** 17 VCF files (Chr00 typically excluded from imputation)

**Note:** BEAGLE mode automatically expects 17 chromosomes, while default mode expects 18.

---

### 2. Extracted Metrics (TSV Columns)

#### Non-BEAGLE Mode: `variant_site_metrics.tsv`

| Column | Description | Source |
|--------|-------------|--------|
| CHROM | Chromosome | VCF CHROM |
| POS | Position | VCF POS |
| QUAL | Variant quality score | VCF QUAL |
| **QD** | Quality by depth | INFO/QD |
| AC | Allele count | INFO/AC |
| AF | Allele frequency | INFO/AF |
| INBREEDING_COEFF | Inbreeding coefficient | INFO/InbreedingCoeff |
| EXCESS_HET | Excess heterozygosity | INFO/ExcessHet |
| **MQ** | Mapping quality | INFO/MQ |
| **MEAN_DEPTH** | Mean read depth across samples | Calculated from FORMAT/DP |
| CALL_RATE | Proportion of samples with called genotypes | Calculated |
| MISSING_RATE | Proportion of missing genotypes | Calculated |
| HETEROZYGOUS_RATE | Proportion of heterozygous calls | Calculated |
| DP_NON_MISSING | Number of samples with depth data | Calculated |
| CALLED_GENOTYPES | Count of called genotypes | Calculated |
| MISSING_GENOTYPES | Count of missing genotypes | Calculated |
| TOTAL_GENOTYPES | Total genotypes (samples × ploidy) | Calculated |
| HETEROZYGOUS_COUNT | Count of heterozygous genotypes | Calculated |

**Total columns:** 18

#### BEAGLE Mode: `variant_site_metrics.tsv`

| Column | Description | Source |
|--------|-------------|--------|
| CHROM | Chromosome | VCF CHROM |
| POS | Position | VCF POS |
| QUAL | Variant quality score | VCF QUAL |
| **AF** | Allele frequency (from imputation) | INFO/AF |
| **DR2** | Dosage R² (imputation quality) | INFO/DR2 |
| **IMP** | Imputation flag (1=imputed, 0=genotyped) | INFO/IMP |
| CALL_RATE | Proportion of samples with called genotypes | Calculated |
| MISSING_RATE | Proportion of missing genotypes | Calculated |
| HETEROZYGOUS_RATE | Proportion of heterozygous calls | Calculated |
| CALLED_GENOTYPES | Count of called genotypes | Calculated |
| MISSING_GENOTYPES | Count of missing genotypes | Calculated |
| TOTAL_GENOTYPES | Total genotypes (samples × ploidy) | Calculated |
| HETEROZYGOUS_COUNT | Count of heterozygous genotypes | Calculated |

**Total columns:** 13

**Key Differences:**
- ❌ No **MEAN_DEPTH, QD, MQ** (not meaningful for imputed data)
- ✅ Added **DR2** (imputation quality metric)
- ✅ Added **IMP** (distinguishes imputed vs genotyped sites)

---

### 3. Generated Plots

#### Non-BEAGLE Mode (Full QC Suite)

| Plot Type | File Pattern | Purpose | Depends On |
|-----------|--------------|---------|------------|
| **Depth vs Position** | `depth_pdfs/Chr*.pdf` | Read depth distribution | DP |
| **Missingness vs Position** | `missingness_plots/Chr*_missingness_vs_position.png` | Missing data patterns | Genotypes |
| **Depth vs Missingness** | `depth_vs_missingness/Chr*_depth_vs_missingness.png` | Depth-missingness correlation | DP |
| **Combined Missingness** | `all_chromosomes_missingness_vs_depth.png` | Overall missing data | All |
| **Site Quality** | Site quality distribution plots | QUAL score distribution | QUAL |
| **Heterozygosity** | Heterozygosity plots | Hardy-Weinberg assessment | Genotypes |
| **Quality by Depth** | QD distribution plots | QD metric assessment | QD |
| **Allele Frequency** | AF distribution plots | MAF spectrum | AF |
| **Call Rate Heatmap** | Call rate heatmaps | Sample-level completeness | Genotypes |

**Total:** ~14 plot types (depth-dependent + general QC)

#### BEAGLE Mode (Imputation QC Suite)

| Plot Type | File Pattern | Purpose | Depends On |
|-----------|--------------|---------|------------|
| **Missingness vs Position** | `missingness_plots/Chr*_missingness_vs_position.png` | Post-imputation missing data | Genotypes |
| **Combined Missingness** | `all_chromosomes_missingness_vs_depth.png` | Overall missing data | All |
| **Site Quality** | Site quality distribution plots | QUAL score distribution | QUAL |
| **Heterozygosity** | Heterozygosity plots | Hardy-Weinberg assessment | Genotypes |
| **Allele Frequency** | AF distribution plots | Imputed AF spectrum | AF |
| **Call Rate Heatmap** | Call rate heatmaps | Sample-level completeness | Genotypes |

**Total:** ~8 plot types (general QC only)

**Plots NOT Generated in BEAGLE Mode:**
- ❌ Depth vs Position (no DP in imputed VCFs)
- ❌ Depth vs Missingness (no DP in imputed VCFs)
- ❌ Quality by Depth (QD not meaningful for imputed data)

**Why skip depth plots?**
- Imputed genotypes don't have read depth information
- DR2 (dosage R²) is the quality metric for imputation, not depth

---

### 4. Usage Examples

#### Non-BEAGLE Mode (Default)

```bash
# Interactive wrapper
bash wrappers/interactive/step1d_interactive.sh --qc

# SLURM submission
bash wrappers/sbatch/step1d_submit.sh Apple_2024 /data/gatk_vcfs --qc

# Direct call
export VCF_DIR="/data/gatk_vcfs"
bash modules/step1d/templates/master_vcf_analysis.sh --qc
```

**Expected Input:**
```
/data/gatk_vcfs/
├── Chr00.vcf.gz
├── Chr01.vcf.gz
├── ...
└── Chr17.vcf.gz
```

#### BEAGLE Mode

```bash
# Interactive wrapper (add --beagle flag)
bash wrappers/interactive/step1d_interactive.sh --beagle --qc

# SLURM submission
bash wrappers/sbatch/step1d_submit.sh Apple_2024 /data/beagle_vcfs --beagle --qc

# Direct call
export VCF_DIR="/data/beagle_vcfs"
bash modules/step1d/templates/master_vcf_analysis.sh --beagle --qc
```

**Expected Input:**
```
/data/beagle_vcfs/
├── Chr01.vcf.gz  # No Chr00
├── Chr02.vcf.gz
├── ...
└── Chr17.vcf.gz
```

---

### 5. Log Messages

#### Non-BEAGLE Mode
```
[step1d] ℹ VCF Directory: /data/gatk_vcfs
[step1d] ℹ Selected mode -> QC: true, PCA: false, Duplicate-check: false
[step1d] ℹ Expected chromosomes: Chr00 Chr01 Chr02 ... Chr17 (18 files)
[step1d] ℹ Extracting site metrics with depth (DP), quality (QD), and mapping quality (MQ)...
```

#### BEAGLE Mode
```
[step1d] ℹ VCF Directory: /data/beagle_vcfs
[step1d] ℹ Beagle mode enabled: using AF/DR2/IMP metrics and skipping depth-dependent outputs.
[step1d] ℹ Selected mode -> QC: true, PCA: false, Duplicate-check: false
[step1d] ℹ Expected chromosomes: Chr01 Chr02 ... Chr17 (17 files, Chr00 excluded)
[step1d] ℹ Extracting site metrics with imputation quality (DR2) and allele frequency (AF)...
[step1d] ℹ Skipping depth plots (mean depth metrics unavailable in BEAGLE mode)
```

---

### 6. Metric Extraction Logic

#### Non-BEAGLE: bcftools query

```bash
bcftools query \
  -f '%CHROM\t%POS\t%QUAL\t%INFO/QD\t%INFO/AC\t%INFO/AF\t%INFO/InbreedingCoeff\t%INFO/ExcessHet\t%INFO/MQ\t[%GT:%DP\t]\n' \
  input.vcf.gz
```

**Extracts:**
- Variant-level: QUAL, QD, AC, AF, InbreedingCoeff, ExcessHet, MQ
- Sample-level: GT (genotype), DP (depth)

#### BEAGLE: bcftools query

```bash
bcftools query \
  -f '%CHROM\t%POS\t%QUAL\t%INFO/AF\t%INFO/DR2\t%INFO/IMP\t[%GT\t]\n' \
  input.vcf.gz
```

**Extracts:**
- Variant-level: QUAL, AF, DR2, IMP
- Sample-level: GT (genotype only, no DP)

---

### 7. INFO Tags Required

#### Non-BEAGLE VCFs (GATK Output)

**Required INFO tags:**
```
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##INFO=<ID=QD,Number=1,Type=Float,Description="Quality by Depth">
##INFO=<ID=MQ,Number=1,Type=Float,Description="Mapping Quality">
##INFO=<ID=InbreedingCoeff,Number=1,Type=Float,Description="Inbreeding coefficient">
##INFO=<ID=ExcessHet,Number=1,Type=Float,Description="Excess heterozygosity">
```

**Required FORMAT tags:**
```
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
```

#### BEAGLE VCFs (Imputed Output)

**Required INFO tags:**
```
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##INFO=<ID=DR2,Number=1,Type=Float,Description="Dosage R-squared (imputation quality)">
##INFO=<ID=IMP,Number=0,Type=Flag,Description="Imputed site">
```

**Required FORMAT tags:**
```
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
```

**Key Difference:** BEAGLE VCFs typically don't have DP (depth) in FORMAT tags because imputation doesn't require read data.

---

### 8. Validation Logic

Step1D automatically validates that existing metrics match the mode:

#### Switching from Non-BEAGLE to BEAGLE

```bash
# First run without --beagle
bash step1d_interactive.sh --qc
# Creates: variant_site_metrics.tsv with MEAN_DEPTH, QD, MQ columns

# Second run WITH --beagle
bash step1d_interactive.sh --beagle --qc
# ⚠️  Detects MEAN_DEPTH column exists
# ⚠️  Regenerates metrics for BEAGLE mode (replaces existing TSV)
```

**Log output:**
```
[step1d] ⚠️  Existing site metrics include depth columns but --beagle mode requested. Regenerating metrics.
```

#### Switching from BEAGLE to Non-BEAGLE

```bash
# First run WITH --beagle
bash step1d_interactive.sh --beagle --qc
# Creates: variant_site_metrics.tsv with AF, DR2, IMP columns

# Second run without --beagle
bash step1d_interactive.sh --qc
# ⚠️  Detects DR2 column exists (BEAGLE-specific)
# ⚠️  Regenerates metrics for non-BEAGLE mode
```

**Log output:**
```
[step1d] ⚠️  Existing site metrics appear to have been generated in --beagle mode. Regenerating metrics for full dataset.
```

**Protection:** Step1D prevents mixing metrics from different modes by auto-detecting and regenerating when mode changes.

---

### 9. PCA Mode Compatibility

Both modes work with PCA:

```bash
# Non-BEAGLE PCA
bash step1d_interactive.sh --PCA

# BEAGLE PCA
bash step1d_interactive.sh --beagle --PCA
```

**No difference in PCA workflow:**
- Both use the same PLINK2 pipeline
- Both apply the same QC filters (geno, mind, maf)
- Both generate the same PCA outputs

**Why?** PCA operates on genotypes (GT field), which both BEAGLE and non-BEAGLE VCFs have.

---

### 10. When to Use Each Mode

#### Use Non-BEAGLE Mode (Default) When:

✅ **Working with raw sequencing data**
- Direct output from GATK HaplotypeCaller
- VCFs from FreeBayes, BCFtools, or other callers
- High-coverage whole genome sequencing (WGS)

✅ **Need depth-based QC**
- Assessing sequencing quality
- Identifying coverage biases
- Detecting low-quality regions

✅ **Evaluating variant calling quality**
- QD (quality by depth) assessment
- Mapping quality (MQ) evaluation
- Variant quality score distribution

**Example scenario:** You just completed GATK variant calling on 150 apple samples and want to assess sequencing and calling quality before downstream analysis.

#### Use BEAGLE Mode When:

✅ **Working with imputed data**
- VCFs processed through BEAGLE software
- Low-coverage sequencing data that was imputed
- Genotype likelihoods converted to hard calls

✅ **Assessing imputation quality**
- DR2 (dosage R²) metric evaluation
- Identifying poorly imputed regions
- Distinguishing imputed vs genotyped sites

✅ **Post-imputation QC**
- Checking imputation completeness
- Validating imputed genotype quality
- Comparing imputed vs original allele frequencies

**Example scenario:** You imputed 50 low-coverage apple samples using a reference panel and want to assess imputation quality before PCA and GWAS.

---

## Practical Examples

### Scenario 1: GATK VCFs (Use Default Mode)

**Your pipeline:**
```bash
# Step 1: Variant calling with GATK
gatk HaplotypeCaller -R reference.fa -I sample1.bam -O sample1.vcf.gz

# Step 2: Joint genotyping
gatk GenotypeGVCFs -R reference.fa -V cohort.vcf.gz -O cohort_genotyped.vcf.gz

# Step 3: Split by chromosome
bcftools view -r Chr01 cohort_genotyped.vcf.gz -Oz -o Chr01.vcf.gz
# ... repeat for all chromosomes ...

# Step 4: QC with Step1D (DEFAULT MODE)
bash step1d_interactive.sh --qc
```

**Why default mode?**
- VCFs have DP (depth) in FORMAT
- VCFs have QD, MQ, etc. in INFO
- No imputation was performed

### Scenario 2: BEAGLE Imputed VCFs (Use BEAGLE Mode)

**Your pipeline:**
```bash
# Step 1: Variant calling with low coverage
gatk HaplotypeCaller -R reference.fa -I sample1.bam -O sample1.vcf.gz --emit-ref-confidence GVCF

# Step 2: Joint genotyping
gatk GenotypeGVCFs -R reference.fa -V cohort.vcf.gz -O cohort_genotyped.vcf.gz

# Step 3: Split and prepare for BEAGLE
for chr in Chr{01..17}; do
    bcftools view -r $chr cohort_genotyped.vcf.gz -Oz -o ${chr}.vcf.gz
done

# Step 4: Run BEAGLE imputation
for chr in Chr{01..17}; do
    java -jar beagle.jar \
        gt=${chr}.vcf.gz \
        ref=reference_panel_${chr}.vcf.gz \
        out=${chr}_imputed
done

# Step 5: QC with Step1D (BEAGLE MODE)
bash step1d_interactive.sh --beagle --qc
```

**Why BEAGLE mode?**
- VCFs have DR2, IMP in INFO (from BEAGLE)
- VCFs don't have meaningful DP (imputed genotypes)
- Want to assess imputation quality (DR2 scores)

---

## Common Questions

### Q: Can I use --beagle with non-imputed VCFs?

**A:** Technically yes, but you'll get errors or missing values:
- ❌ DR2 column will be "NA" (tag doesn't exist)
- ❌ IMP column will be "0" (tag doesn't exist)
- ❌ Plots will be generated but less meaningful

**Recommendation:** Don't use `--beagle` unless your VCFs are actually imputed.

### Q: What if my BEAGLE VCFs still have DP tags?

**A:** Some BEAGLE workflows preserve original DP for genotyped sites:
- ✅ Use `--beagle` mode (focuses on DR2/IMP)
- ℹ️ DP will be ignored in metrics extraction
- ℹ️ Depth plots will be skipped

### Q: Can I compare metrics between modes?

**A:** Not directly - they measure different things:
- Non-BEAGLE: Sequencing/calling quality (depth-based)
- BEAGLE: Imputation quality (DR2-based)

### Q: Will PCA results differ between modes?

**A:** No, if the genotypes (GT) are the same:
- PCA uses GT field only
- Mode affects QC metrics, not genotypes
- PCA results depend on variant filtering, not mode

### Q: What about other imputation tools (Minimac, Impute)?

**A:** BEAGLE mode works if the tool outputs similar tags:
- ✅ Need: AF, DR2 (or similar R² metric), IMP (or imputation flag)
- ⚠️ Tag names must match exactly: `INFO/DR2`, `INFO/IMP`
- 💡 Consider renaming tags: `bcftools annotate --rename-annots`

---

## Summary Table

| Feature | Non-BEAGLE | BEAGLE |
|---------|------------|--------|
| **Typical source** | GATK, FreeBayes, BCFtools | BEAGLE imputation software |
| **Data type** | Raw genotypes | Imputed genotypes |
| **Chromosomes** | 18 (Chr00-Chr17) | 17 (Chr01-Chr17) |
| **Key metrics** | DP, QD, MQ | DR2, IMP |
| **Plots** | 14+ types | 8 types |
| **Depth plots** | ✅ Yes | ❌ No |
| **Imputation QC** | ❌ No | ✅ Yes (DR2) |
| **Use for** | Sequencing QC | Imputation QC |
| **PCA compatible** | ✅ Yes | ✅ Yes |
| **Flag** | (none - default) | `--beagle` |

---

## Related Documentation

- `README_FIRST.md` - Main Step1D documentation
- `FUNCTION_MAP.md` - Complete flag reference
- `OUTPUT_REVIEW.md` - All output files explained
- `USAGE_DECISION_TREE.md` - When to use which mode

---

**Last Updated:** February 3, 2026
