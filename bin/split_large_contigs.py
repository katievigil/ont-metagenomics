#!/usr/bin/env python3
"""Write every contig >= --min-length to its own FASTA file (input for CheckM2)
and record its length and circularity.

Circularity is read from the metaMDBG FASTA header (e.g. 'circular=true');
see the README for this assumption."""
import argparse
import os
import re
import sys

CIRCULAR = re.compile(r"circular\s*[=:]\s*(true|yes|1)\b", re.IGNORECASE)


def safe_stem(name):
    return re.sub(r"[^A-Za-z0-9._-]", "_", name)


def read_fasta(path):
    """Yield (name, header, sequence); name is the first whitespace-delimited token."""
    name, header, chunks = None, None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    yield name, header, "".join(chunks)
                header = line[1:]
                name = header.split()[0]
                chunks = []
            elif line:
                chunks.append(line)
    if name is not None:
        yield name, header, "".join(chunks)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contigs", required=True)
    parser.add_argument("--min-length", type=int, required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--info", required=True)
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    written = 0
    with open(args.info, "w") as info:
        info.write("contig\tfile_stem\tlength\tcircular\n")
        for name, header, seq in read_fasta(args.contigs):
            if len(seq) < args.min_length:
                continue
            stem = safe_stem(name)
            circular = int(bool(CIRCULAR.search(header)))
            with open(os.path.join(args.outdir, stem + ".fa"), "w") as out:
                out.write(f">{name}\n{seq}\n")
            info.write(f"{name}\t{stem}\t{len(seq)}\t{circular}\n")
            written += 1
    print(f"{written} contigs >= {args.min_length} bp written to {args.outdir}", file=sys.stderr)


if __name__ == "__main__":
    main()
