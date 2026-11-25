#!/usr/bin/env python3
"""Per-sample summary tables, mirroring the app note's reported metrics:

  run_summary.tsv  read, assembly, MAG tier and AMR counts (cf. app note Table 1)
  amr_loci.tsv     AMR loci with contig, MAG and chromosome/plasmid localisation
  mag_table.tsv    MAGs with quality tier, GTDB-Tk taxonomy and AMR loci count"""
import argparse
import csv
import glob
import os
import sys
from collections import Counter


def read_fasta_lengths(path):
    lengths, current, started = [], 0, False
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if started:
                    lengths.append(current)
                current, started = 0, True
            else:
                current += len(line.strip())
    if started:
        lengths.append(current)
    return lengths


def n50(lengths):
    total, acc = sum(lengths), 0
    for length in sorted(lengths, reverse=True):
        acc += length
        if acc * 2 >= total:
            return length
    return 0


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path, header, rows):
    with open(path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for row in rows:
            fh.write("\t".join(str(x) for x in row) + "\n")


def load_taxonomy(gtdbtk_dir):
    taxonomy = {}
    for path in glob.glob(os.path.join(gtdbtk_dir, "*summary.tsv")):
        for row in read_tsv(path):
            taxonomy[row["user_genome"]] = row["classification"]
    return taxonomy


def species_from_gtdb(classification):
    for rank in classification.split(";"):
        rank = rank.strip()
        if rank.startswith("s__") and len(rank) > 3:
            return rank[3:]
    return "unclassified"


def classify_viralverify(prediction):
    p = prediction.lower()
    if p.startswith("non") or "chrom" in p:
        return "chromosomal"
    if "plasmid" in p:
        return "plasmid"
    return "uncertain"  # viral and uncertain calls are grouped as uncertain (app note)


def load_localisation(viralverify_dir):
    """Contig -> chromosomal / plasmid / uncertain, from the viralVerify prediction table."""
    for path in sorted(glob.glob(os.path.join(viralverify_dir, "*"))):
        if os.path.isdir(path):
            continue
        with open(path) as fh:
            lines = [line.rstrip("\n") for line in fh if line.strip()]
        if not lines or "prediction" not in lines[0].lower():
            continue
        header = lines[0]
        delim = "," if header.count(",") > header.count("\t") else "\t"
        cols = [c.strip().strip('"') for c in header.split(delim)]
        pred_idx = next(i for i, c in enumerate(cols) if "prediction" in c.lower())
        calls = {}
        for line in lines[1:]:
            fields = [f.strip().strip('"') for f in line.split(delim)]
            if len(fields) > pred_idx:
                calls[fields[0].split()[0]] = classify_viralverify(fields[pred_idx])
        return calls
    print("WARNING: no viralVerify prediction table found; plasmid calls unavailable",
          file=sys.stderr)
    return {}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--read-stats", required=True)
    parser.add_argument("--contigs", required=True)
    parser.add_argument("--amr", required=True)
    parser.add_argument("--viralverify", required=True)
    parser.add_argument("--mag-summary", required=True)
    parser.add_argument("--contig2mag", required=True)
    parser.add_argument("--gtdbtk", required=True)
    args = parser.parse_args()

    stats = read_tsv(args.read_stats)
    reads = stats[0] if stats else {}
    lengths = read_fasta_lengths(args.contigs)
    contig2mag = {row["contig"]: row["mag_id"] for row in read_tsv(args.contig2mag)}
    mags = read_tsv(args.mag_summary)
    taxonomy = load_taxonomy(args.gtdbtk)
    localisation = load_localisation(args.viralverify)
    amr_hits = read_tsv(args.amr)

    # AMR loci
    location_counts, mag_amr_counts = Counter(), Counter()
    amr_contigs, amr_rows, n_fragmented = set(), [], 0
    for hit in amr_hits:
        contig = hit.get("Contig id", "")
        location = localisation.get(contig, "unknown")
        mag_id = contig2mag.get(contig, "unbinned")
        method = hit.get("Method", "")
        status = "fragmented" if "PARTIAL_CONTIG_END" in method else "complete"
        location_counts[location] += 1
        mag_amr_counts[mag_id] += 1
        amr_contigs.add(contig)
        n_fragmented += int(status == "fragmented")
        amr_rows.append([args.sample, contig, mag_id, hit.get("Element symbol", ""),
                         hit.get("Element name", ""), hit.get("Class", ""),
                         hit.get("Subclass", ""), method, status, location])
    write_tsv("amr_loci.tsv",
              ["sample", "contig", "mag_id", "element_symbol", "element_name", "class",
               "subclass", "method", "status", "location"],
              amr_rows)

    # MAG table
    tiers, species, mag_rows = Counter(), set(), []
    for mag in mags:
        mag_id = mag["mag_id"]
        classification = taxonomy.get(mag_id, "unclassified")
        sp = species_from_gtdb(classification)
        if sp != "unclassified":
            species.add(sp)
        tiers[mag["tier"]] += 1
        mag_rows.append([mag_id, mag["tier"], mag["source"], mag["completeness"],
                         mag["contamination"], mag["length_bp"], classification, sp,
                         mag_amr_counts.get(mag_id, 0)])
    write_tsv("mag_table.tsv",
              ["mag_id", "tier", "source", "completeness", "contamination", "length_bp",
               "gtdb_classification", "species", "n_amr_loci"],
              mag_rows)

    # Run summary
    closed, sc_hq = tiers.get("closed", 0), tiers.get("SC_HQ", 0)
    hq, mq = tiers.get("HQ", 0), tiers.get("MQ", 0)
    summary = [
        ("reads_after_filter", reads.get("num_seqs", "NA")),
        ("bases_after_filter", reads.get("sum_len", "NA")),
        ("read_N50", reads.get("N50", "NA")),
        ("mean_read_quality", reads.get("AvgQual", "NA")),
        ("contigs", len(lengths)),
        ("assembly_bp", sum(lengths)),
        ("contig_N50", n50(lengths)),
        ("contigs_ge_1Mb", sum(1 for length in lengths if length >= 1_000_000)),
        ("closed_MAGs", closed),
        ("SC_HQ_MAGs_linear", sc_hq),
        ("SC_HQ_MAGs_total", closed + sc_hq),
        ("HQ_MAGs_total", closed + sc_hq + hq),
        ("MQ_plus_HQ_MAGs", closed + sc_hq + hq + mq),
        ("species_classified", len(species)),
        ("AMR_loci", len(amr_rows)),
        ("AMR_loci_chromosomal", location_counts["chromosomal"]),
        ("AMR_loci_plasmid", location_counts["plasmid"]),
        ("AMR_loci_uncertain_or_unknown", location_counts["uncertain"] + location_counts["unknown"]),
        ("AMR_loci_fragmented", n_fragmented),
        ("AMR_contigs", len(amr_contigs)),
    ]
    write_tsv("run_summary.tsv", ["metric", "value"], summary)
    print(f"{args.sample}: {len(mags)} MAGs, {len(amr_rows)} AMR loci", file=sys.stderr)


if __name__ == "__main__":
    main()
