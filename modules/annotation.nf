// Taxonomy (GTDB-Tk), AMR detection (AMRFinderPlus), plasmid calls (viralVerify), per-sample summary.

process GTDBTK_CLASSIFY {
    tag "${sample}"
    label 'process_high'
    conda 'bioconda::gtdbtk=2.4.0'
    publishDir "${params.outdir}/${sample}/taxonomy", mode: 'copy'

    input:
    tuple val(sample), path(mags)

    output:
    tuple val(sample), path('gtdbtk_out'), emit: taxonomy

    script:
    """
    export GTDBTK_DATA_PATH="${params.gtdbtk_db}"

    if [ -n "\$(ls -A ${mags})" ]; then
        gtdbtk classify_wf \\
            --genome_dir ${mags} \\
            --extension fa \\
            --out_dir gtdbtk_out \\
            --cpus ${task.cpus}
    else
        mkdir -p gtdbtk_out
        echo "no MAGs to classify" > gtdbtk_out/no_mags.txt
    fi
    """
}

process AMRFINDER {
    tag "${sample}"
    label 'process_medium'
    conda 'bioconda::ncbi-amrfinderplus=4.0.3'
    publishDir "${params.outdir}/${sample}/amr", mode: 'copy'

    input:
    tuple val(sample), path(contigs)

    output:
    tuple val(sample), path('amr_contigs.tsv'), emit: amr

    script:
    def db = params.amrfinder_db ? "--database \"${params.amrfinder_db}\"" : ''
    """
    awk '/^>/ { print \$1; next } { print }' "${contigs}" > contigs.clean.fa

    amrfinder \\
        --nucleotide contigs.clean.fa \\
        --output amr_contigs.tsv \\
        --threads ${task.cpus} \\
        ${db}
    """
}

process VIRALVERIFY {
    tag "${sample}"
    label 'process_medium'
    conda 'bioconda::viralverify=1.1'
    publishDir "${params.outdir}/${sample}/plasmids", mode: 'copy'

    input:
    tuple val(sample), path(contigs)

    output:
    tuple val(sample), path('viralverify'), emit: results

    script:
    """
    awk '/^>/ { print \$1; next } { print }' "${contigs}" > contigs.clean.fa

    viralverify \\
        -f contigs.clean.fa \\
        -o viralverify \\
        -hmm "${params.viralverify_hmm}" \\
        -p
    """
}

process SUMMARISE {
    tag "${sample}"
    label 'process_low'
    conda 'conda-forge::python=3.11'
    publishDir "${params.outdir}/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(read_stats), path(contigs), path(amr), path(viralverify), path(mag_summary), path(contig2mag), path(gtdbtk)

    output:
    path 'run_summary.tsv'
    path 'amr_loci.tsv'
    path 'mag_table.tsv'

    script:
    """
    summarise_run.py \\
        --sample "${sample}" \\
        --read-stats "${read_stats}" \\
        --contigs "${contigs}" \\
        --amr "${amr}" \\
        --viralverify "${viralverify}" \\
        --mag-summary "${mag_summary}" \\
        --contig2mag "${contig2mag}" \\
        --gtdbtk "${gtdbtk}"
    """
}
