#!/usr/bin/env python3

from argparse import ArgumentParser
import glob
import json
import os
import sys
from typing import List

from Bio import SeqIO


def main() -> None:
    parser = ArgumentParser()
    parser.add_argument("files", nargs="*",
                        help="file(s) to check, defaults to all files in retired when empty")
    parser.add_argument("--genbank-dir", "-g", default="genbanks",
                        help="Directory containing the GenBank files (default: %(default)s).")
    args = parser.parse_args()

    if not args.files:
        args.files = sorted(glob.glob(os.path.join("retired", "*.json")))

    run(args.files, args.genbank_dir)


def run(files: List[str], genbank_dir: str) -> None:
    for file in files:
        with open(file, 'r', encoding="utf-8") as handle:
            data = json.load(handle)
        accession: str = data["cluster"]["loci"]["accession"]
        if accession.startswith("MIBIG.BGC"):
            continue

        filename = os.path.join(genbank_dir, f"{accession}.gbk")
        record = SeqIO.read(filename, "genbank")
        cds_count = gene_count = 0
        for feature in record.features:
            if feature.type == "CDS":
                cds_count += 1
            elif feature.type == "gene":
                gene_count += 1

        if cds_count == 1 or gene_count == 1:
            mibig_acc = data['cluster']['mibig_accession']
            print(f"{mibig_acc}: Only {cds_count} CDSes/ {gene_count} genes", file=sys.stderr)


if __name__ == "__main__":
    main()
