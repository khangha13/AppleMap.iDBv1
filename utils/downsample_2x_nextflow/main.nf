#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.accessions = params.accessions ?: null
params.bam_dir = params.bam_dir ?: null
params.outdir = params.outdir ?: null
params.bam_suffix = params.bam_suffix ?: '.bam'
params.target_depth = params.target_depth ?: 2
params.seed = params.seed ?: 67
params.chrom_regex = params.chrom_regex ?: '^Chr(0[1-9]|1[0-7])$'

if (!params.accessions) {
    error 'Missing required parameter: --accessions'
}

if (!params.bam_dir) {
    error 'Missing required parameter: --bam_dir'
}

if (!params.outdir) {
    error 'Missing required parameter: --outdir'
}

def accessions_file = file(params.accessions)
def bam_dir = file(params.bam_dir)

if (!accessions_file.exists()) {
    error "Accession list not found: ${params.accessions}"
}

if (!bam_dir.exists()) {
    error "BAM directory not found: ${params.bam_dir}"
}

workflow {
    Channel
        .fromPath(accessions_file)
        .splitText()
        .map { line -> line.trim() }
        .filter { accession -> accession }
        .map { accession ->
            def bam_path = file("${bam_dir}/${accession}${params.bam_suffix}")
            if (!bam_path.exists()) {
                error "BAM file not found for accession ${accession}: ${bam_path}"
            }
            tuple(accession, bam_path)
        }
        | DOWNSAMPLE_BAM
}

process DOWNSAMPLE_BAM {
    tag "${accession}"

    cpus 4
    memory '8 GB'
    time '24h'
    module 'samtools/1.18-gcc-12.3.0'

    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(accession), path(bam)

    output:
    tuple val(accession), path("${accession}_2x.bam"), path("${accession}_2x.bam.bai")

    script:
    """
    set -euo pipefail

    TARGET_DEPTH="${params.target_depth}"
    SEED="${params.seed}"
    CHROM_REGEX='${params.chrom_regex}'
    OUTPUT_BAM="${accession}_2x.bam"

    echo "=========================================="
    echo "Processing accession: ${accession}"
    echo "Input BAM: ${bam}"
    echo "Output BAM: \${OUTPUT_BAM}"
    echo "Target depth: \${TARGET_DEPTH}x"
    echo "Threads: ${task.cpus}"
    echo "=========================================="

    if [ ! -f "${bam}" ]; then
        echo "ERROR: Staged BAM not found: ${bam}"
        exit 1
    fi

    if [ ! -f "${bam}.bai" ]; then
        echo "BAM index not found in task work directory. Creating index: ${bam}.bai"
        samtools index -@ ${task.cpus} "${bam}"
    fi

    echo "Estimating length-weighted mean depth over Chr01-Chr17..."
    depth=\$(
        samtools coverage "${bam}" | awk -v chr_re="\${CHROM_REGEX}" '
            \$1 ~ chr_re {
                len = \$3 - \$2 + 1
                depth_sum += \$7 * len
                length_sum += len
            }
            END {
                if (length_sum > 0) {
                    printf "%.6f", depth_sum / length_sum
                } else {
                    print "NA"
                }
            }'
    )

    if [[ ! "\${depth}" =~ ^[0-9]+([.][0-9]+)?\$ ]]; then
        echo "ERROR: Could not calculate depth over Chr01-Chr17 for ${accession} (depth=\${depth})"
        exit 1
    fi

    echo "Observed depth across Chr01-Chr17: \${depth}x"

    if awk -v depth="\${depth}" -v target="\${TARGET_DEPTH}" 'BEGIN { exit !(depth <= target) }'; then
        echo "ERROR: ${accession} is already at or below \${TARGET_DEPTH}x (observed \${depth}x). Not downsampling."
        exit 1
    fi

    fraction_arg=\$(awk -v seed="\${SEED}" -v target="\${TARGET_DEPTH}" -v depth="\${depth}" 'BEGIN { printf "%.6f", seed + target / depth }')
    keep_fraction=\$(awk -v target="\${TARGET_DEPTH}" -v depth="\${depth}" 'BEGIN { printf "%.6f", target / depth }')

    echo "Downsampling keep fraction: \${keep_fraction}"
    echo "samtools view -s argument: \${fraction_arg}"

    samtools view -@ ${task.cpus} -b -s "\${fraction_arg}" "${bam}" -o "\${OUTPUT_BAM}"
    samtools index -@ ${task.cpus} "\${OUTPUT_BAM}"

    echo "Downsampling complete for ${accession}"
    echo "Output BAM: \${OUTPUT_BAM}"
    echo "Output index: \${OUTPUT_BAM}.bai"
    """
}
