# ONT metagenomics pipeline (Nextflow)

This pipeline implements the workflow of Oxford Nanopore application note AN_1282 (v2), *Oxford Nanopore sequencing provides superior MAG recovery and strain-level resolution from a complex microbiome* (ZymoBIOMICS Fecal Reference). It runs from reads to metagenome-assembled genomes (MAGs), taxonomy, AMR loci, and chromosome/plasmid localisation of AMR genes.

## Workflow

| Step | Process(es) | Tool (version in the note) | Notes |
|---|---|---|---|
| Basecalling (POD5 folders only) | `BASECALL_DORADO` | Dorado; default SUP v5.0.0 + 6mA + 4mC_5mC, as in the note | Choose the model from the dropdown. Needs a GPU and a Dorado container or binary |
| Read filtering | `FILTER_READS` | Samtools 1.18 | Mean read Q ≥ 10 and length ≥ 1 kb |
| Merge runs | `CONCAT_READS` | SeqKit | All reads for a sample are combined before assembly |
| Assembly | `METAMDBG_ASSEMBLE` | metaMDBG 1.0 | `--in-ont` |
| Read mapping | `MAP_READS` | minimap2 2.27, Samtools | `lr:hq`, secondary alignments excluded |
| SC HQ MAGs | `SPLIT_LARGE_CONTIGS`, `CHECKM2` (large), `PARTITION_SC_HQ` | CheckM2 1.0.2 | Contigs ≥ 1 Mb. Linear ≥ 95 % or circular ≥ 90 % completeness |
| Binning | `SEMIBIN2`, `METABAT2` | SemiBin2 2.1.0 (human_gut, long read), MetaBAT2 2.17 | Run on contigs that are not SC HQ MAGs |
| Consolidation | `DAS_TOOL` | DAS Tool 1.1.7 | |
| HQ/MQ MAGs | `CHECKM2` (bins), `BUILD_MAG_SET` | CheckM2 | HQ: ≥ 90 % / ≤ 5 %. MQ: ≥ 50 % / ≤ 10 % |
| Taxonomy | `GTDBTK_CLASSIFY` | GTDB-Tk 2.4.0, GTDB R220 | |
| AMR genes | `AMRFINDER` | AMRFinderPlus 4.0.3 | Fragmented if the method is `PARTIAL_CONTIG_END*` |
| Plasmids | `VIRALVERIFY` | viralVerify 1.1 (`-p`) | Viral and uncertain calls are grouped as uncertain |
| Summary | `SUMMARISE` | Python | Run metrics, AMR loci table, MAG table |

The note's comparisons with Illumina and PacBio HiFi data, its gggenomes synteny plot, and the Strainy/devider haplotype analyses are not included.

## Requirements

- Nextflow 23.10 or later
- One of: Conda/Mamba (`-profile conda`), Docker, Singularity, or Wave (`-profile wave`)
- A GPU for Dorado, only if you start from POD5 folders
- Large memory: the note's data (about 200 Gb per PromethION flow cell) needs hundreds of GB of RAM for metaMDBG. Test on a subsample first.
- Databases:
  - CheckM2: `checkm2 database --download --path <dir>`
  - GTDB-Tk: GTDB R220 reference data (about 110 GB, unpacked)
  - AMRFinderPlus: `amrfinder_update --force_update --database <dir>`
  - viralVerify: the `nbc_hmms.hmm` file

## Usage

```bash
nextflow run katievigil/ont-metagenomics -r main -profile conda \
    --fastq reads.fastq.gz \
    --checkm2_db /db/checkm2 \
    --gtdbtk_db /db/release220 \
    --viralverify_hmm /db/nbc_hmms.hmm \
    --amrfinder_db /db/amrfinder \
    --outdir results
```

Use `-resume` to restart after a failure.

## Inputs

| Parameter | Description |
|---|---|
| `--fastq` | One of: a single FASTQ(.gz) or BAM file; a folder of FASTQ or BAM files, which are combined into one sample; or a folder of POD5 files, which are basecalled first |
| `--sample` | Optional sample name for the output folder. Defaults to the file or folder name |

To analyse two flow cells as one sample, put both files in one folder and pass the folder to `--fastq`.

## Outputs

```
results/<sample>/
  reads/            <sample>.read_stats.tsv
  assembly/         contigs.fa (metaMDBG)
  mags/             mags/*.fa, mag_summary.tsv, contig2mag.tsv, sc_hq_summary.tsv
  taxonomy/         gtdbtk_out/  (GTDB-Tk summary tables)
  amr/              amr_contigs.tsv  (AMRFinderPlus)
  plasmids/         viralverify/     (viralVerify output)
  run_summary.tsv   read, assembly, MAG tier and AMR counts
  amr_loci.tsv      AMR loci: contig, MAG, element, class, chromosome/plasmid
  mag_table.tsv     MAG tier, completeness, contamination, GTDB taxonomy, AMR count
results/pipeline_info/  timeline, report and trace
```

## Points to check before a production run

The pipeline was written from the note's Methods and from the tools' documentation. It has not been run end to end on real data. Check these against your installed versions:

1. **Samtools filter syntax.** `FILTER_READS` uses `mean(qual)` and `qlen` in `samtools view -e`.
2. **Dorado models.** Choose the model with `--basecall_model` and the modification models with `--basecall_mods`. The dropdown lists the DNA models from the Oxford Nanopore model list. Run `dorado download --list` in the Dorado container to confirm the installed version supports them, and confirm that `model,mods` is accepted with `dorado basecaller --help`.
3. **SemiBin2 output.** `SEMIBIN2` expects `contig_bins.tsv` in the output folder.
4. **Circular contigs.** `split_large_contigs.py` treats a contig as circular when its metaMDBG header contains `circular=true`. Closed MAGs depend on this flag.
5. **viralVerify output.** `summarise_run.py` finds the table with a `Prediction` column in the viralVerify output folder.
6. **Containers.** The default `conda` profile installs bioconda packages with the note's versions where bioconda has them.
7. **Resources.** Adjust the labels in `nextflow.config` for your machine.

## Running in EPI2ME Desktop

1. Push this repository to GitHub as `katievigil/ont-metagenomics`.
2. In EPI2ME Desktop, go to Workflows, click Import workflow, and enter `https://github.com/katievigil/ont-metagenomics`.
3. Docker must be installed and running.

EPI2ME reads `nextflow.config` and `nextflow_schema.json` to build the launch form, so keep both in sync when you add or rename parameters.
