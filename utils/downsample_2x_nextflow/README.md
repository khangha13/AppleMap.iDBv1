# Nextflow 2x BAM Downsampling

Downsample a flat directory of BAM files to 2x coverage using one Nextflow task per accession. The workflow also writes the complementary reads that were not selected into a separate remainder BAM.

This workflow can run either locally or on Bunya through the `bunya` SLURM profile.

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

Existing BAM indexes can be either BAI or CSI:

```text
<accession><bam_suffix>.bai
<accession><bam_suffix>.csi
```

If neither index exists, the workflow creates an index in the task work directory. It tries BAI first and falls back to CSI if BAI creation fails.

## Run Locally

From the repository root:

```bash
nextflow run utils/downsample_2x_nextflow/main.nf \
  --accessions accessions.txt \
  --bam_dir /path/to/bams \
  --bam_suffix .recal.bam \
  --outdir /path/to/downsampled_2x \
  -resume
```

Local mode does not use SLURM. Nextflow runs each task in the current session using the default local executor. Use this for small tests, development, or an interactive compute session where `samtools` is already available.

If your BAMs are named `<accession>.bam`, omit `--bam_suffix .recal.bam`.

## Run On Bunya With SLURM

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

The `-profile bunya` option tells Nextflow to use the `bunya` profile in `nextflow.config`. That profile submits each accession as a separate SLURM job with these settings:

- Executor: `slurm`
- Account: `a_qaafi_cas`
- Partition: `general`
- Resources per accession: `4 CPUs`, `8 GB`, `24h`
- Module: `samtools/1.18-gcc-12.3.0`

Use this mode for production runs on Bunya.

`--outdir` is required. Outputs are not published beside the source BAMs unless you explicitly set `--outdir` to that location.

## Parameters

| Parameter | Required | Default | Description |
| --- | --- | --- | --- |
| `--accessions` | Yes | None | Text file with one accession per line. |
| `--bam_dir` | Yes | None | Flat directory containing `<accession><bam_suffix>`. |
| `--outdir` | Yes | None | Directory where downsampled BAMs, remainder BAMs, and indexes are published. |
| `--bam_suffix` | No | `.bam` | BAM filename suffix appended to each accession. Use `.recal.bam` for files like `ID.recal.bam`. |
| `--target_depth` | No | `2` | Target downsampling depth. |
| `--seed` | No | `67` | Seed used for `samtools view -s`. |
| `--chrom_regex` | No | `^Chr(0[1-9]\|1[0-7])$` | Chromosomes used for depth estimation. |

## Output

For each accession, the workflow publishes the 2x downsampled BAM and the complementary remainder BAM:

```text
<accession>_2x.bam
<accession>_2x.bam.bai or <accession>_2x.bam.csi
<accession>_remainder_<actualCoverage>x.bam
<accession>_remainder_<actualCoverage>x.bam.bai or <accession>_remainder_<actualCoverage>x.bam.csi
<accession>.downsampling_metrics.tsv
```

The workflow also publishes a combined report:

```text
downsampling_metrics.tsv
```

The remainder BAM contains reads from the input BAM that were not selected by `samtools view -s` for the 2x BAM. The two BAM outputs are complementary splits of the source BAM.

Example:

```text
SRR000001_2x.bam
SRR000001_2x.bam.bai
SRR000001_remainder_6.238000x.bam
SRR000001_remainder_6.238000x.bam.bai
SRR000001.downsampling_metrics.tsv
```

Metrics columns:

```text
sample  original_depth  target_depth  keep_fraction  downsampled_depth  remainder_depth  downsampled_bam  remainder_bam
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

After the split, the workflow calculates the actual coverage of the remainder BAM using the same Chr01-Chr17 method. That measured value is written into the remainder filename. It is usually close to `original depth - target depth`, but can differ slightly because random subsampling is stochastic.

## Failure Conditions

A task fails if:

- The required accession list, BAM directory, or output directory parameter is missing.
- `<accession><bam_suffix>` is missing from `--bam_dir`.
- Mean depth cannot be calculated over `Chr01` to `Chr17`.
- Observed depth is already at or below the target depth.
- Remainder BAM coverage cannot be calculated after splitting.

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
