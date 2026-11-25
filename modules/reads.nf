// Read handling: basecalling (POD5 only), quality/length filtering, merging runs per sample.

process BASECALL_DORADO {
    tag "${sample}/${run}"
    label 'gpu'
    container params.dorado_container

    input:
    tuple val(sample), val(run), path(pod5_dir)

    output:
    tuple val(sample), val(run), path("${sample}.${run}.basecalls.bam"), emit: bam

    script:
    def model = params.basecall_mods == 'none' ? params.basecall_model : "${params.basecall_model},${params.basecall_mods}"
    """
    dorado basecaller "${model}" "${pod5_dir}" \\
        --device "${params.dorado_device}" \\
        --recursive \\
        > "${sample}.${run}.basecalls.bam"
    """
}

process FILTER_READS {
    tag "${sample}/${run}"
    label 'process_low'
    conda 'bioconda::samtools=1.18'

    input:
    tuple val(sample), val(run), path(reads)

    output:
    tuple val(sample), path("${sample}.${run}.filtered.fastq.gz")

    script:
    """
    set -euo pipefail

    # Dorado BAM is used directly; FASTQ is imported as an unaligned BAM first
    case "${reads}" in
        *.bam) IN_BAM="${reads}" ;;
        *)     samtools import -0 "${reads}" -o imported.bam
               IN_BAM=imported.bam ;;
    esac

    # App note: mean read quality >= 10 and read length >= 1 kb
    samtools view -b \\
        -e 'qlen >= ${params.min_len} && mean(qual) >= ${params.min_qual}' \\
        -o filtered.bam "\${IN_BAM}"

    samtools fastq filtered.bam | gzip -c > "${sample}.${run}.filtered.fastq.gz"
    """
}

process CONCAT_READS {
    tag "${sample}"
    label 'process_low'
    conda 'bioconda::seqkit>=2.0'
    publishDir "${params.outdir}/${sample}/reads", mode: 'copy', pattern: '*.read_stats.tsv'

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("${sample}.filtered.fastq.gz"), emit: reads
    tuple val(sample), path("${sample}.read_stats.tsv"), emit: stats

    script:
    """
    cat ${reads} > "${sample}.filtered.fastq.gz"
    seqkit stats -T -a "${sample}.filtered.fastq.gz" > "${sample}.read_stats.tsv"
    """
}
