#!/usr/bin/env python3
"""Assemble the final MAG set (app note downstream analysis):

  * single-contig MAGs from contigs >= 1 Mb (tier 'closed' or 'SC_HQ'),
    as produced by partition_sc_hq.py
  * binned MAGs from DAS Tool consolidated bins that reach HQ or MQ quality
    in CheckM2 (HQ: completeness >= 90 and contamination <= 5;
                MQ: completeness >= 50 and contamination <= 10)

Writes one FASTA per MAG, a MAG summary table and a contig-to-MAG table."""
import argparse
import csv
import os
import shutil
import sys

FASTA_EXTENSIONS = (".fa", ".fasta", ".fna")


def read_fasta(path):
    name, chunks = None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            elif line:
                chunks.append(line)
    if name is not None:
        yield name, "".join(chunks)


def load_checkm2(path):
    with open(path) as fh:
        return {row["Name"]: (float(row["Completeness"]), float(row["Contamination"]))
                for row in csv.DictReader(fh, delimiter="\t")}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sc-hq-dir", required=True)
    parser.add_argument("--sc-hq-summary", required=True)
    parser.add_argument("--bins-dir", required=True)
    parser.add_argument("--bin-checkm2", required=True)
    parser.add_argument("--hq-completeness", type=float, default=90.0)
    parser.add_argument("--hq-contamination", type=float, default=5.0)
    parser.add_argument("--mq-completeness", type=float, default=50.0)
    parser.add_argument("--mq-contamination", type=float, default=10.0)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--mag-summary", required=True)
    parser.add_argument("--contig2mag", required=True)
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    summary_rows, contig_map = [], []

    # 1. Single-contig MAGs
    with open(args.sc_hq_summary) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            mag_id = "SCHQ_" + row["file_stem"]
            shutil.copy(os.path.join(args.sc_hq_dir, row["file_stem"] + ".fa"),
                        os.path.join(args.outdir, mag_id + ".fa"))
            summary_rows.append({
                "mag_id": mag_id,
                "tier": row["tier"],
                "source": "single_contig",
                "completeness": row["completeness"],
                "contamination": row["contamination"],
                "length_bp": row["length"],
                "n_contigs": 1,
            })
            contig_map.append((row["contig"], mag_id))

    # 2. Binned MAGs passing HQ or MQ thresholds
    quality = load_checkm2(args.bin_checkm2)
    for fn in sorted(os.listdir(args.bins_dir)):
        if not fn.endswith(FASTA_EXTENSIONS):
            continue
        stem = os.path.splitext(fn)[0]
        if stem not in quality:
            print(f"WARNING: no CheckM2 result for {fn}; skipping", file=sys.stderr)
            continue
        completeness, contamination = quality[stem]
        if completeness >= args.hq_completeness and contamination <= args.hq_contamination:
            tier = "HQ"
        elif completeness >= args.mq_completeness and contamination <= args.mq_contamination:
            tier = "MQ"
        else:
            continue

        mag_id = "BIN_" + stem
        seqs = list(read_fasta(os.path.join(args.bins_dir, fn)))
        with open(os.path.join(args.outdir, mag_id + ".fa"), "w") as out:
            for name, seq in seqs:
                out.write(f">{name}\n{seq}\n")
                contig_map.append((name, mag_id))
        summary_rows.append({
            "mag_id": mag_id,
            "tier": tier,
            "source": "binned",
            "completeness": completeness,
            "contamination": contamination,
            "length_bp": sum(len(seq) for _, seq in seqs),
            "n_contigs": len(seqs),
        })

    fields = ["mag_id", "tier", "source", "completeness", "contamination", "length_bp", "n_contigs"]
    with open(args.mag_summary, "w") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(summary_rows)

    with open(args.contig2mag, "w") as fh:
        fh.write("contig\tmag_id\n")
        for contig, mag_id in contig_map:
            fh.write(f"{contig}\t{mag_id}\n")

    counts = {}
    for row in summary_rows:
        counts[row["tier"]] = counts.get(row["tier"], 0) + 1
    tiers = ", ".join(f"{tier}={n}" for tier, n in sorted(counts.items())) or "none"
    print(f"MAG set: {len(summary_rows)} MAGs ({tiers})", file=sys.stderr)


if __name__ == "__main__":
    main()
