// Quality assessment, binning (SemiBin2 + MetaBAT2), DAS Tool consolidation and final MAG set.

process CHECKM2 {
    tag "${sample}"
    label 'process_medium'
    conda 'bioconda::checkm2=1.0.2'

    input:
    tuple val(sample), path(genomes)

    output:
    tuple val(sample), path('checkm2/quality_report.tsv'), emit: report

    script:
    """
    if [ -n "\$(ls -A ${genomes})" ]; then
        checkm2 predict \\
            --threads ${task.cpus} \\
            --input ${genomes} \\
            --output-directory checkm2 \\
            --database_path "${params.checkm2_db}" \\
            -x fa
    else
        mkdir -p checkm2
        printf 'Name\\tCompleteness\\tContamination\\n' > checkm2/quality_report.tsv
    fi
    """
}

process SEMIBIN2 {
    tag "${sample}"
    label 'process_high'
    conda 'bioconda::semibin=2.1.0'

    input:
    tuple val(sample), path(bam), path(bai), path(contigs)

    output:
    tuple val(sample), path('semibin_contig2bin.tsv'), emit: contig2bin

    script:
    """
    SemiBin2 single_easy_bin \\
        -i "${contigs}" \\
        -b "${bam}" \\
        -o semibin \\
        --environment ${params.semibin_env} \\
        --sequencing-type long_read \\
        --min-len ${params.semibin_min_len} \\
        --threads ${task.cpus}

    awk -F'\\t' 'NF >= 2 && \$1 != "contig_name" && \$1 != "contig" { print \$1 "\\t" \$2 }' \\
        semibin/contig_bins.tsv > semibin_contig2bin.tsv
    """
}

process METABAT2 {
    tag "${sample}"
    label 'process_medium'
    conda 'bioconda::metabat2=2.17 bioconda::das_tool=1.1.7'

    input:
    tuple val(sample), path(bam), path(bai), path(contigs)

    output:
    tuple val(sample), path('metabat2_contig2bin.tsv'), emit: contig2bin

    script:
    """
    jgi_summarize_bam_contig_depths --outputDepth depth.txt "${bam}"

    mkdir -p metabat2
    metabat2 \\
        -i "${contigs}" \\
        -a depth.txt \\
        -o metabat2/bin \\
        -m ${params.metabat_min_len} \\
        -t ${task.cpus}

    Fasta_to_Contig2Bin.sh -i metabat2 -e fa > metabat2_contig2bin.tsv
    """
}

process DAS_TOOL {
    tag "${sample}"
    label 'process_medium'
    conda 'bioconda::das_tool=1.1.7'

    input:
    tuple val(sample), path(semibin), path(metabat2), path(contigs)

    output:
    tuple val(sample), path('dastool_DASTool_bins'), emit: bins

    script:
    """
    DAS_Tool \\
        -i "${semibin},${metabat2}" \\
        -l SemiBin2,MetaBAT2 \\
        -c "${contigs}" \\
        -o dastool \\
        --write_bins 1 \\
        --threads ${task.cpus}

    mkdir -p dastool_DASTool_bins
    """
}

process BUILD_MAG_SET {
    tag "${sample}"
    label 'process_low'
    conda 'conda-forge::python=3.11'
    publishDir "${params.outdir}/${sample}/mags", mode: 'copy'

    input:
    tuple val(sample), path(sc_hq_dir), path(sc_hq_summary), path(bins_dir), path(checkm2_report)

    output:
    tuple val(sample), path('mags'), emit: mags
    tuple val(sample), path('mag_summary.tsv'), path('contig2mag.tsv'), emit: tables

    script:
    """
    build_mag_set.py \\
        --sc-hq-dir "${sc_hq_dir}" \\
        --sc-hq-summary "${sc_hq_summary}" \\
        --bins-dir "${bins_dir}" \\
        --bin-checkm2 "${checkm2_report}" \\
        --hq-completeness ${params.hq_completeness} \\
        --hq-contamination ${params.hq_contamination} \\
        --mq-completeness ${params.mq_completeness} \\
        --mq-contamination ${params.mq_contamination} \\
        --outdir mags \\
        --mag-summary mag_summary.tsv \\
        --contig2mag contig2mag.tsv
    """
}
