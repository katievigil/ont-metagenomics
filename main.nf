#!/usr/bin/env nextflow
/*
 * ONT metagenomics pipeline
 *
 * Follows the workflow of Oxford Nanopore application note AN_1282 (v2),
 * "Oxford Nanopore sequencing provides superior MAG recovery and strain-level
 * resolution from a complex microbiome" (ZymoBIOMICS Fecal Reference):
 *
 *   POD5 -> Dorado basecalling                       (optional, per run)
 *        -> Samtools quality/length filter           (per run)
 *        -> concatenate runs per sample
 *        -> metaMDBG assembly -> minimap2 read mapping
 *        -> single-contig HQ MAGs from contigs >= 1 Mb (CheckM2)
 *        -> SemiBin2 + MetaBAT2 on remaining contigs -> DAS Tool -> CheckM2
 *        -> final MAG set (SC HQ + HQ/MQ bins)
 *        -> GTDB-Tk taxonomy, AMRFinderPlus, viralVerify (chromosome vs plasmid)
 */

nextflow.enable.dsl = 2

include { BASECALL_DORADO; FILTER_READS; CONCAT_READS } from './modules/reads'
include { METAMDBG_ASSEMBLE; MAP_READS; SPLIT_LARGE_CONTIGS; PARTITION_SC_HQ } from './modules/assembly'
include { CHECKM2 as CHECKM2_LARGE; CHECKM2 as CHECKM2_BINS; SEMIBIN2; METABAT2; DAS_TOOL; BUILD_MAG_SET } from './modules/binning'
include { GTDBTK_CLASSIFY; AMRFINDER; VIRALVERIFY; SUMMARISE } from './modules/annotation'

def helpMessage() {
    log.info """
    ONT metagenomics pipeline (after Oxford Nanopore application note AN_1282)

    Usage:
      nextflow run main.nf -profile conda --samplesheet assets/samplesheet.csv \\
          --checkm2_db /db/checkm2 --gtdbtk_db /db/gtdbtk_r220 --viralverify_hmm /db/nbc_hmms.hmm

    Required:
      --samplesheet        CSV with columns sample,run,reads (see README)
      --checkm2_db         CheckM2 database directory
      --gtdbtk_db          GTDB-Tk reference data directory (GTDB R220)
      --viralverify_hmm    viralVerify HMM file (nbc_hmms.hmm)

    Optional:
      --outdir             output directory (default: results)
      --amrfinder_db       AMRFinderPlus database directory
      --large_contig_size  contigs at or above this size are assessed individually (default: 1000000)
    """.stripIndent()
}

workflow {

    if (params.help) {
        helpMessage()
        exit 0
    }
    if (!params.samplesheet)     exit 1, 'ERROR: --samplesheet is required (see --help)'
    if (!params.checkm2_db)      exit 1, 'ERROR: --checkm2_db is required'
    if (!params.gtdbtk_db)       exit 1, 'ERROR: --gtdbtk_db is required'
    if (!params.viralverify_hmm) exit 1, 'ERROR: --viralverify_hmm is required'

    // Inputs: one row per sequencing run; rows sharing a sample id are combined
    ch_runs = Channel
        .fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample, row.run, file(row.reads, checkIfExists: true)) }

    ch_runs_by_type = ch_runs.branch { sample, run, reads ->
        pod5:  reads.isDirectory()
        other: true
    }

    // Per run: basecall POD5 directories (if given), then filter reads
    BASECALL_DORADO(ch_runs_by_type.pod5)
    FILTER_READS(BASECALL_DORADO.out.bam.mix(ch_runs_by_type.other))

    // Per sample: combine runs
    CONCAT_READS(FILTER_READS.out.groupTuple())
    ch_reads = CONCAT_READS.out.reads

    // Assembly and read mapping
    METAMDBG_ASSEMBLE(ch_reads)
    ch_contigs = METAMDBG_ASSEMBLE.out.contigs
    MAP_READS(ch_contigs.join(ch_reads))

    // Single-contig HQ MAGs from contigs >= 1 Mb
    SPLIT_LARGE_CONTIGS(ch_contigs)
    CHECKM2_LARGE(SPLIT_LARGE_CONTIGS.out.large)
    PARTITION_SC_HQ(
        ch_contigs
            .join(SPLIT_LARGE_CONTIGS.out.info)
            .join(CHECKM2_LARGE.out.report)
    )

    // Binning of every contig that is not a SC HQ MAG
    ch_bin_in = MAP_READS.out.join(PARTITION_SC_HQ.out.remaining)
    SEMIBIN2(ch_bin_in)
    METABAT2(ch_bin_in)
    DAS_TOOL(
        SEMIBIN2.out.contig2bin
            .join(METABAT2.out.contig2bin)
            .join(PARTITION_SC_HQ.out.remaining)
    )
    CHECKM2_BINS(DAS_TOOL.out.bins)

    // Final MAG set: SC HQ MAGs plus HQ/MQ bins
    BUILD_MAG_SET(
        PARTITION_SC_HQ.out.sc_hq
            .join(DAS_TOOL.out.bins)
            .join(CHECKM2_BINS.out.report)
    )

    // Annotation
    GTDBTK_CLASSIFY(BUILD_MAG_SET.out.mags)
    AMRFINDER(ch_contigs)
    VIRALVERIFY(ch_contigs)

    SUMMARISE(
        CONCAT_READS.out.stats
            .join(ch_contigs)
            .join(AMRFINDER.out.amr)
            .join(VIRALVERIFY.out.results)
            .join(BUILD_MAG_SET.out.tables)
            .join(GTDBTK_CLASSIFY.out.taxonomy)
    )
}
