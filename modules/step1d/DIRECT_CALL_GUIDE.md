# Direct Script Call with Custom VCF Filename and Directory

## Quick Command

```bash
# Set your custom VCF directory and filename
export VCF_DIR="/path/to/your/vcf/directory"
export VCF_INCLUDE_FILENAMES="your_file.vcf.gz"
export WORK_DIR="${VCF_DIR}"

# Run PCA-only mode directly
bash /path/to/GATK_Pipeline_KH_v1/modules/step1d/templates/master_vcf_analysis.sh --PCA
```

## Complete Example

Here's a complete example with all the necessary settings:

```bash
#!/bin/bash

# 1. Set pipeline root (adjust to your path)
PIPELINE_ROOT="/Users/khangha/Library/CloudStorage/OneDrive-TheUniversityofQueensland/Code/Apple_GATK_2025/KH/GATK_Pipeline_KH_v1"

# 2. Set your custom VCF directory
export VCF_DIR="/custom/path/to/your/vcfs"

# 3. Set your VCF filename(s)
# Single file:
export VCF_INCLUDE_FILENAMES="my_custom_file.vcf.gz"

# Multiple files (space-separated):
# export VCF_INCLUDE_FILENAMES="file1.vcf.gz file2.vcf.gz file3.vcf.gz"

# 4. Set working directory (usually same as VCF_DIR)
export WORK_DIR="${VCF_DIR}"

# 5. Optional: Set R scripts directory (auto-detected if not set)
# export R_SCRIPTS_DIR="${PIPELINE_ROOT}/modules/step1d/Rscripts"

# 6. Run the script
bash "${PIPELINE_ROOT}/modules/step1d/templates/master_vcf_analysis.sh" --PCA
```

## Single-Line Command (Copy-Paste Ready)

Replace the paths and filename with your actual values:

```bash
export VCF_DIR="/your/custom/path" && export VCF_INCLUDE_FILENAMES="your_file.vcf.gz" && export WORK_DIR="${VCF_DIR}" && bash /path/to/GATK_Pipeline_KH_v1/modules/step1d/templates/master_vcf_analysis.sh --PCA
```

## Multiple VCF Files

If you have multiple VCF files in the same directory:

```bash
export VCF_DIR="/path/to/vcfs"
export VCF_INCLUDE_FILENAMES="file1.vcf.gz file2.vcf.gz file3.vcf.gz"
export WORK_DIR="${VCF_DIR}"

bash /path/to/GATK_Pipeline_KH_v1/modules/step1d/templates/master_vcf_analysis.sh --PCA
```

## Important Notes

1. **VCF_INCLUDE_FILENAMES**: 
   - Use **space-separated** filenames (not commas)
   - Include the full filename with extension (e.g., `file.vcf.gz` or `file.vcf`)
   - Files must exist in `VCF_DIR`

2. **VCF_DIR**: 
   - Use **absolute path** (full path starting with `/`)
   - Directory must exist and contain your VCF file(s)

3. **WORK_DIR**: 
   - Usually same as `VCF_DIR`
   - This is where output files will be created

4. **Output Location**: 
   - Results will be in: `${WORK_DIR}/pca_analysis/`
   - Default directory name is `pca_analysis` (can be changed with `STEP1D_PCA_DIR`)

## Verification Before Running

Check that your file exists:

```bash
# Set your values
export VCF_DIR="/path/to/your/vcfs"
export VCF_INCLUDE_FILENAMES="your_file.vcf.gz"

# Verify file exists
ls -lh "${VCF_DIR}/${VCF_INCLUDE_FILENAMES}"
```

## Common Issues

### Issue: "No VCF files specified for analysis"
**Solution:** Make sure `VCF_INCLUDE_FILENAMES` is set correctly:
```bash
# Check if variable is set
echo "VCF_INCLUDE_FILENAMES=${VCF_INCLUDE_FILENAMES}"

# Make sure filename matches exactly (case-sensitive)
ls -la "${VCF_DIR}/${VCF_INCLUDE_FILENAMES}"
```

### Issue: "VCF directory not found"
**Solution:** Use absolute path and verify it exists:
```bash
# Use absolute path (starting with /)
export VCF_DIR="/absolute/path/to/vcfs"

# Verify directory exists
test -d "${VCF_DIR}" && echo "Directory exists" || echo "Directory NOT found"
```

### Issue: File not found
**Solution:** Check filename spelling and extension:
```bash
# List files in directory to see exact names
ls -la "${VCF_DIR}"

# Make sure VCF_INCLUDE_FILENAMES matches exactly
```

## Optional Customizations

### Change PCA output directory name:
```bash
export STEP1D_PCA_DIR="my_custom_pca_folder"
```

### Remove relatives before PCA:
```bash
export STEP1D_REMOVE_RELATIVES=true
# Note: Requires --pca-only mode (already enabled)
```

### For Beagle-imputed VCFs:
```bash
bash "${MASTER_SCRIPT}" --beagle --pca-only
```

## Output Files

After running, you'll find:
```
${WORK_DIR}/pca_analysis/
├── pca.eigenvec          # Sample coordinates
├── pca.eigenval          # Eigenvalues
├── pca_PC1_PC2.png       # PCA plot
├── pca_PC1_PC3.png       # PCA plot
├── pca_PC2_PC3.png       # PCA plot
└── pca_scree.png         # Scree plot
```

## Complete Working Example

```bash
#!/bin/bash
# Complete example - replace paths with your actual values

# Pipeline location
PIPELINE_ROOT="/Users/khangha/Library/CloudStorage/OneDrive-TheUniversityofQueensland/Code/Apple_GATK_2025/KH/GATK_Pipeline_KH_v1"

# Your VCF location
export VCF_DIR="/scratch/user/myproject/custom_vcfs"
export VCF_INCLUDE_FILENAMES="my_genome_data.vcf.gz"
export WORK_DIR="${VCF_DIR}"
export STEP1D_PCA_ONLY=true

# Verify file exists
if [ ! -f "${VCF_DIR}/${VCF_INCLUDE_FILENAMES}" ]; then
    echo "Error: VCF file not found: ${VCF_DIR}/${VCF_INCLUDE_FILENAMES}"
    exit 1
fi

# Run PCA
bash "${PIPELINE_ROOT}/modules/step1d/templates/master_vcf_analysis.sh" --pca-only

echo "Done! Check results in: ${WORK_DIR}/pca_analysis/"
```
