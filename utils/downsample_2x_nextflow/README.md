# Nextflow 2x BAM Downsampling

Downsample a flat directory of BAM files to 2x coverage using one Nextflow task per accession.

This workflow is intended for Bunya and uses:

- SLURM executor
- `a_qaafi_cas` account
- `general` partition
- `samtools/1.18-gcc-12.3.0`

## Inputs

The accession list must contain one accession per line:

```text
SRR000001
SRR000002
SRR000003
```

The BAM directory must be flat and contain one BAM per accession. By default, the workflow expects `<accession>.bam`:

```text
/path/to/bams/
|-- SRR000001.bam
|-- SRR000002.bam
`-- SRR000003.bam
```

For BAMs named like `<accession>.recal.bam`, pass `--bam_suffix .recal.bam`.

## Run On Bunya

From the repository root:

```bash
module load nextflow

nextflow run utils/downsample_2x_nextflow/main.nf \
  -profile bunya \
  --accessions accessions.txt \
  --bam_dir /path/to/bams \
  --bam_suffix .recal.bam \
  --outdir /path/to/downsampled_2x \
  -resume
```

`--outdir` is required. Outputs are not published beside the source BAMs unless you explicitly set `--outdir` to that location.

## Parameters

| Parameter | Required | Default | Description |
| --- | --- | --- | --- |
| `--accessions` | Yes | None | Text file with one accession per line. |
| `--bam_dir` | Yes | None | Flat directory containing `<accession><bam_suffix>`. |
| `--outdir` | Yes | None | Directory where downsampled BAMs and indexes are published. |
| `--bam_suffix` | No | `.bam` | BAM filename suffix appended to each accession. Use `.recal.bam` for files like `ID.recal.bam`. |
| `--target_depth` | No | `2` | Target downsampling depth. |
| `--seed` | No | `67` | Seed used for `samtools view -s`. |
| `--chrom_regex` | No | `^Chr(0[1-9]\|1[0-7])$` | Chromosomes used for depth estimation. |

## Output

For each accession, the workflow publishes:

```text
<accession>_2x.bam
<accession>_2x.bam.bai
```

## Depth Calculation

Depth is calculated with `samtools coverage` and a length-weighted mean across `Chr01` to `Chr17`.

The keep fraction is:

```text
target_depth / observed_depth
```

The `samtools view -s` argument is:

```text
seed + target_depth / observed_depth
```

For example, with seed `67`, target `2x`, and observed depth `8x`, the `samtools view -s` argument is `67.250000`.

## Failure Conditions

A task fails if:

- The required accession list, BAM directory, or output directory parameter is missing.
- `<accession><bam_suffix>` is missing from `--bam_dir`.
- Mean depth cannot be calculated over `Chr01` to `Chr17`.
- Observed depth is already at or below the target depth.

These failures are intentional so low-coverage BAMs are flagged instead of silently copied or overwritten.

## Resume

Use `-resume` to reuse completed Nextflow work directories after fixing inputs or rerunning the workflow:

```bash
nextflow run utils/downsample_2x_nextflow/main.nf \
  -profile bunya \
  --accessions accessions.txt \
  --bam_dir /path/to/bams \
  --bam_suffix .recal.bam \
  --outdir /path/to/downsampled_2x \
  -resume
```
