#!/usr/bin/env python3
"""Identify single-contig high-quality (SC HQ) MAGs among contigs >= 1 Mb and
write every other contig to --remaining for binning.

Criteria (app note): linear contigs with CheckM2 completeness >= 95 %,
or circular contigs with completeness >= 90 %. Circular SC HQ MAGs are
tiered 'closed'; linear ones are tiered 'SC_HQ'."""
import argparse
import csv
import os
import sys


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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contigs", required=True)
    parser.add_argument("--contig-info", required=True)
    parser.add_argument("--checkm2", required=True)
    parser.add_argument("--linear-completeness", type=float, default=95.0)
    parser.add_argument("--circular-completeness", type=float, default=90.0)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--remaining", required=True)
    args = parser.parse_args()

    with open(args.checkm2) as fh:
        quality = {row["Name"]: (float(row["Completeness"]), float(row["Contamination"]))
                   for row in csv.DictReader(fh, delimiter="\t")}

    with open(args.contig_info) as fh:
        info = {row["contig"]: row for row in csv.DictReader(fh, delimiter="\t")}

    sc_hq = {}
    for contig, row in info.items():
        if row["file_stem"] not in quality:
            continue
        completeness, contamination = quality[row["file_stem"]]
        circular = row["circular"] == "1"
        threshold = args.circular_completeness if circular else args.linear_completeness
        if completeness >= threshold:
            sc_hq[contig] = {
                "contig": contig,
                "file_stem": row["file_stem"],
                "length": row["length"],
                "circular": row["circular"],
                "completeness": completeness,
                "contamination": contamination,
                "tier": "closed" if circular else "SC_HQ",
            }

    os.makedirs(args.outdir, exist_ok=True)
    n_remaining = 0
    with open(args.remaining, "w") as rem:
        for name, seq in read_fasta(args.contigs):
            if name in sc_hq:
                stem = sc_hq[name]["file_stem"]
                with open(os.path.join(args.outdir, stem + ".fa"), "w") as out:
                    out.write(f">{name}\n{seq}\n")
            else:
                rem.write(f">{name}\n{seq}\n")
                n_remaining += 1

    fields = ["contig", "file_stem", "length", "circular", "completeness", "contamination", "tier"]
    with open(args.summary, "w") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(sc_hq.values())

    closed = sum(1 for r in sc_hq.values() if r["tier"] == "closed")
    print(f"SC HQ MAGs: {len(sc_hq)} ({closed} closed); contigs left for binning: {n_remaining}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
