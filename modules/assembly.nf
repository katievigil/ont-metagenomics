// Assembly, read mapping, and single-contig HQ MAG identification.

process METAMDBG_ASSEMBLE {
    tag "${sample}"
    label 'process_high'
    conda 'bioconda::metamdbg=1.0'
    publishDir "${params.outdir}/${sample}/assembly", mode: 'copy', pattern: 'contigs.fa'

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path('contigs.fa'), emit: contigs

    script:
    """
    metaMDBG asm \\
        --out-dir metamdbg \\
        --in-ont "${reads}" \\
        --threads ${task.cpus}

    gunzip -c metamdbg/contigs.fasta.gz > contigs.fa
    """
}

process MAP_READS {
    tag "${sample}"
    label 'process_medium'
    conda 'bioconda::minimap2=2.27 bioconda::samtools=1.18'

    input:
    tuple val(sample), path(contigs), path(reads)

    output:
    tuple val(sample), path("${sample}.mapped.sorted.bam"), path("${sample}.mapped.sorted.bam.bai")

    script:
    """
    minimap2 -t ${task.cpus} -ax lr:hq --secondary=no "${contigs}" "${reads}" \\
        | samtools sort -@ ${task.cpus} -o "${sample}.mapped.sorted.bam" -
    samtools index "${sample}.mapped.sorted.bam"
    """
}

process SPLIT_LARGE_CONTIGS {
    tag "${sample}"
    label 'process_low'
    conda 'conda-forge::python=3.11'

    input:
    tuple val(sample), path(contigs)

    output:
    tuple val(sample), path('large_contigs'), emit: large
    tuple val(sample), path('contig_info.tsv'), emit: info

    script:
    """
    split_large_contigs.py \\
        --contigs "${contigs}" \\
        --min-length ${params.large_contig_size} \\
        --outdir large_contigs \\
        --info contig_info.tsv
    """
}

process PARTITION_SC_HQ {
    tag "${sample}"
    label 'process_low'
    conda 'conda-forge::python=3.11'
    publishDir "${params.outdir}/${sample}/mags", mode: 'copy', pattern: 'sc_hq_summary.tsv'

    input:
    tuple val(sample), path(contigs), path(contig_info), path(checkm2_report)

    output:
    tuple val(sample), path('sc_hq'), path('sc_hq_summary.tsv'), emit: sc_hq
    tuple val(sample), path('remaining.fa'), emit: remaining

    script:
    """
    partition_sc_hq.py \\
        --contigs "${contigs}" \\
        --contig-info "${contig_info}" \\
        --checkm2 "${checkm2_report}" \\
        --linear-completeness ${params.sc_hq_linear_completeness} \\
        --circular-completeness ${params.sc_hq_circular_completeness} \\
        --outdir sc_hq \\
        --summary sc_hq_summary.tsv \\
        --remaining remaining.fa
    """
}
