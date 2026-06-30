#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.accessions = params.accessions ?: null
params.bam_dir = params.bam_dir ?: null
params.outdir = params.outdir ?: null
params.bam_suffix = params.bam_suffix ?: '.bam'
params.target_depth = params.target_depth ?: 2
params.seed = params.seed ?: 67
params.chrom_regex = params.chrom_regex ?: '^Chr(0[1-9]|1[0-7])$'

workflow {
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
            def bai_path = file("${bam_path}.bai")
            def csi_path = file("${bam_path}.csi")
            def index_path = bai_path.exists() ? bai_path.toString() : (csi_path.exists() ? csi_path.toString() : '')
            tuple(accession, bam_path, index_path)
        }
        | DOWNSAMPLE_BAM

    DOWNSAMPLE_BAM.out.metrics.collect() | MERGE_DOWNSAMPLING_METRICS
}

process DOWNSAMPLE_BAM {
    tag "${accession}"

    cpus 4
    memory '8 GB'
    time '24h'
    module 'samtools/1.18-gcc-12.3.0'

    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(accession), path(bam), val(source_index)

    output:
    tuple val(accession),
        path("${accession}_2x.bam"),
        path("${accession}_2x.bam.*"),
        path("${accession}_remainder_*x.bam"),
        path("${accession}_remainder_*x.bam.*")
    path("${accession}.downsampling_metrics.tsv"), emit: metrics

    script:
    """
    set -euo pipefail

    TARGET_DEPTH="${params.target_depth}"
    SEED="${params.seed}"
    CHROM_REGEX='${params.chrom_regex}'
    OUTPUT_BAM="${accession}_2x.bam"
    REMAINDER_TMP_BAM="${accession}_remainder.tmp.bam"
    SOURCE_INDEX='${source_index}'

    calculate_depth() {
        local depth_bam="\$1"
        samtools coverage "\${depth_bam}" | awk -v chr_re="\${CHROM_REGEX}" '
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
    }

    ensure_bam_index() {
        local index_bam="\$1"
        if [ -f "\${index_bam}.bai" ] || [ -f "\${index_bam}.csi" ]; then
            echo "BAM index found for \${index_bam}"
            return 0
        fi

        echo "BAM index not found for \${index_bam}. Creating BAI index."
        if samtools index -@ ${task.cpus} "\${index_bam}"; then
            return 0
        fi

        echo "BAI index creation failed for \${index_bam}. Creating CSI index instead."
        samtools index -@ ${task.cpus} -c "\${index_bam}"
    }

    echo "=========================================="
    echo "Processing accession: ${accession}"
    echo "Input BAM: ${bam}"
    echo "Output BAM: \${OUTPUT_BAM}"
    echo "Remainder temporary BAM: \${REMAINDER_TMP_BAM}"
    echo "Target depth: \${TARGET_DEPTH}x"
    echo "Threads: ${task.cpus}"
    echo "=========================================="

    if [ ! -f "${bam}" ]; then
        echo "ERROR: Staged BAM not found: ${bam}"
        exit 1
    fi

    if [ -n "\${SOURCE_INDEX}" ] && [ -f "\${SOURCE_INDEX}" ]; then
        case "\${SOURCE_INDEX}" in
            *.bai)
                ln -sf "\${SOURCE_INDEX}" "${bam}.bai"
                ;;
            *.csi)
                ln -sf "\${SOURCE_INDEX}" "${bam}.csi"
                ;;
            *)
                echo "WARNING: Unrecognised BAM index extension, ignoring: \${SOURCE_INDEX}"
                ;;
        esac
    fi

    ensure_bam_index "${bam}"

    echo "Estimating length-weighted mean depth over Chr01-Chr17..."
    depth=\$(calculate_depth "${bam}")

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

    samtools view -@ ${task.cpus} -b -s "\${fraction_arg}" -U "\${REMAINDER_TMP_BAM}" -o "\${OUTPUT_BAM}" "${bam}"
    ensure_bam_index "\${OUTPUT_BAM}"

    echo "Estimating actual coverage of downsampled BAM over Chr01-Chr17..."
    downsampled_depth=\$(calculate_depth "\${OUTPUT_BAM}")

    if [[ ! "\${downsampled_depth}" =~ ^[0-9]+([.][0-9]+)?\$ ]]; then
        echo "ERROR: Could not calculate downsampled depth for ${accession} (depth=\${downsampled_depth})"
        exit 1
    fi

    echo "Estimating actual coverage of remainder BAM over Chr01-Chr17..."
    remainder_depth=\$(calculate_depth "\${REMAINDER_TMP_BAM}")

    if [[ ! "\${remainder_depth}" =~ ^[0-9]+([.][0-9]+)?\$ ]]; then
        echo "ERROR: Could not calculate remainder depth for ${accession} (depth=\${remainder_depth})"
        exit 1
    fi

    REMAINDER_BAM="${accession}_remainder_\${remainder_depth}x.bam"
    mv "\${REMAINDER_TMP_BAM}" "\${REMAINDER_BAM}"
    ensure_bam_index "\${REMAINDER_BAM}"

    {
        printf 'sample\\toriginal_depth\\ttarget_depth\\tkeep_fraction\\tdownsampled_depth\\tremainder_depth\\tdownsampled_bam\\tremainder_bam\\n'
        printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' \
            '${accession}' \
            "\${depth}" \
            "\${TARGET_DEPTH}" \
            "\${keep_fraction}" \
            "\${downsampled_depth}" \
            "\${remainder_depth}" \
            "\${OUTPUT_BAM}" \
            "\${REMAINDER_BAM}"
    } > "${accession}.downsampling_metrics.tsv"

    echo "Downsampling complete for ${accession}"
    echo "Output BAM: \${OUTPUT_BAM}"
    echo "Output index: \$(ls \${OUTPUT_BAM}.bai \${OUTPUT_BAM}.csi 2>/dev/null | tr '\\n' ' ')"
    echo "Downsampled observed depth: \${downsampled_depth}x"
    echo "Remainder observed depth: \${remainder_depth}x"
    echo "Remainder BAM: \${REMAINDER_BAM}"
    echo "Remainder index: \$(ls \${REMAINDER_BAM}.bai \${REMAINDER_BAM}.csi 2>/dev/null | tr '\\n' ' ')"
    """
}

process MERGE_DOWNSAMPLING_METRICS {
    tag 'downsampling_metrics'

    publishDir params.outdir, mode: 'copy'

    input:
    path metrics_files

    output:
    path 'downsampling_metrics.tsv'

    script:
    """
    set -euo pipefail

    first_file=\$(ls *.downsampling_metrics.tsv | sort | head -n 1)
    awk 'FNR == 1 && NR != 1 { next } { print }' \$(ls *.downsampling_metrics.tsv | sort) > downsampling_metrics.tsv

    if [ ! -s downsampling_metrics.tsv ]; then
        echo "ERROR: downsampling_metrics.tsv was not created" >&2
        exit 1
    fi

    echo "Merged metrics using header from \${first_file}"
    """
}
