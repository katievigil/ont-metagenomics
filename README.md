# ONT metagenomics pipeline (Nextflow)

This pipeline implements the workflow of Oxford Nanopore application note AN_1282 (v2), *Oxford Nanopore sequencing provides superior MAG recovery and strain-level resolution from a complex microbiome* (ZymoBIOMICS Fecal Reference). It runs from raw reads to metagenome-assembled genomes (MAGs), taxonomy, AMR loci and their chromosome/plasmid localisation.

## Workflow

| Step | Process(es) | Tool (version in the note) | Notes |
|---|---|---|---|
| Basecalling (POD5 inputs only) | `BASECALL_DORADO` | Dorado, SUP v5.0.0 + 6mA + 4mC_5mC | Needs a GPU and a Dorado container or binary |
| Read filtering | `FILTER_READS` | Samtools 1.18 | Mean read Q ≥ 10 and length ≥ 1 kb |
| Merge runs | `CONCAT_READS` | SeqKit | Runs sharing a `sample` id are combined, as the note does for its two flow cells |
| Assembly | `METAMDBG_ASSEMBLE` | metaMDBG 1.0 | `--in-ont` |
| Read mapping | `MAP_READS` | minimap2 2.27, Samtools | `lr:hq`, secondary alignments excluded |
| SC HQ MAGs | `SPLIT_LARGE_CONTIGS`, `CHECKM2` (large), `PARTITION_SC_HQ` | CheckM2 1.0.2 | Contigs ≥ 1 Mb. Linear ≥ 95 % or circular ≥ 90 % completeness |
| Binning | `SEMIBIN2`, `METABAT2` | SemiBin2 2.1.0 (human_gut, long read), MetaBAT2 2.17 | Run on the contigs that are not SC HQ MAGs |
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
- A GPU for Dorado, only if you start from POD5
- Large memory: the note's data (~200 Gb per PromethION flow cell) needs hundreds of GB of RAM for metaMDBG. Test on a subsample first.
- Databases:
  - CheckM2: `checkm2 database --download --path <dir>`
  - GTDB-Tk: GTDB R220 reference data (about 100 GB or more, unpacked)
  - AMRFinderPlus: `amrfinder -u` or `amrfinder_update --force_update --database <dir>`
  - viralVerify: the `nbc_hmms.hmm` file

## Usage

```bash
nextflow run main.nf -profile conda \
    --samplesheet assets/samplesheet.csv \
    --checkm2_db /db/checkm2 \
    --gtdbtk_db /db/gtdbtk_r220 \
    --viralverify_hmm /db/viralverify/nbc_hmms.hmm \
    --amrfinder_db /db/amrfinder \
    --outdir results
```

Use `-resume` to restart after a failure.

## Samplesheet

One row per sequencing run. Rows with the same `sample` are combined before assembly.

| Column | Description |
|---|---|
| `sample` | Sample ID used for output folders |
| `run` | Run ID (for example a flow cell). Must be unique within a sample |
| `reads` | A POD5 directory (basecalled by the pipeline), an unaligned BAM from Dorado, or a FASTQ(.gz) file |

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

The pipeline was written from the note's Methods and from the tools' documentation. It has not been run end to end. Check these against your installed versions:

1. **Samtools filter syntax.** `FILTER_READS` uses `mean(qual)` and `qlen` in `samtools view -e`. Confirm with `samtools view --help` if the filter fails.
2. **Dorado model and modification flags.** The default `sup@v5.0.0,6mA,4mC_5mC` is set with `--dorado_model`. Confirm with `dorado basecaller --help`. The Dorado container tag defaults to `latest`. Pin it.
3. **SemiBin2 output.** `SEMIBIN2` expects `contig_bins.tsv` in the output folder. Adjust the awk line in `modules/binning.nf` if your version writes a different file.
4. **Circular contigs.** `split_large_contigs.py` treats a contig as circular when its metaMDBG header contains `circular=true`. If your metaMDBG version writes this differently, edit `CIRCULAR` in that script. Closed MAGs depend on this flag.
5. **viralVerify output.** `summarise_run.py` finds the table with a `Prediction` column in the viralVerify output folder. Check that the column and contig names match.
6. **Containers.** The default `conda` profile installs bioconda packages with the versions from the note, where bioconda has them. For Docker or Singularity, use `-profile wave`, which builds containers from the conda directives. Without Wave, processes run on the host even when Docker is enabled.
7. **Resources.** Adjust the labels in `nextflow.config` for your cluster.

## Testing

Run on a small subsample first:

1. Subsample the reads (for example `samtools view -s`, or `seqkit sample -p 0.05`).
2. Run with `-profile conda` and a single sample.
3. Check `run_summary.tsv` and `mag_table.tsv`.

## Running in EPI2ME Desktop

1. Push this repository to GitHub as `katievigil/ont-metagenomics`.
2. In EPI2ME Desktop, go to Workflows, click Import workflow, and enter `https://github.com/katievigil/ont-metagenomics`.
3. Docker must be installed and running.

EPI2ME reads `nextflow.config` and `nextflow_schema.json` to build the launch form, so keep both in sync when you add or rename parameters.
# ont-metagenomics
