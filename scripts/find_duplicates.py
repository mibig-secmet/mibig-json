#!/usr/bin/env python3

import argparse
from collections import defaultdict
import glob
import json
import os
import sys
from typing import List, Tuple

END_GUARD = 1e9


def compare(prev: Tuple[int, int, str], current: Tuple[int, int, str], threshold: int) -> bool:
    # exact coordinate match
    if prev[:2] == current[:2]:
        print("entries duplicated:", prev[2], current[2])
        return True

    # if threshold is -1, don't check for overlaps at all
    if threshold < 0:
        return False

    overlap = 0
    if prev[0] <= current[0] <= prev[1]:
        if prev[1] != END_GUARD:
            overlap = prev[1] - current[0]
        elif current[1] == END_GUARD:
            overlap = f"all but {current[0] - prev[0]}"
        else:
            overlap = "almost all"
    elif prev[0] <= current[1] <= prev[1]:
        overlap = current[1] - prev[0]

    if overlap:
        if not isinstance(overlap, int):
            print("entries overlap:", prev[2], current[2], "end positions unknown")
            return True
        if overlap > threshold:
            print("entries overlap:", prev[2], current[2], " by ", overlap, "nucleotides")
            return True

    return False


def run_all(files: List[str], threshold: int) -> bool:
    if not files:
        files = []
        for dir_name in ["data", "pending"]:
            files.extend(glob.glob(os.path.join(dir_name, "*.json")))

    accessions = defaultdict(list)

    for file in sorted(files):
        with open(file) as handle:
            data = json.load(handle)
        loci = data["cluster"]["loci"]
        acc = loci["accession"].rsplit(".", 1)[0]
        start = loci.get("start_coord", 0)
        end = loci.get("end_coord", END_GUARD)

        accessions[acc].append((start, end, data["cluster"]["mibig_accession"]))

    found = False
    for locations in accessions.values():
        locations = sorted(locations)
        prev = locations[0]
        for location in locations[1:]:
            found = compare(prev, location, threshold) or found
            prev = location
    return found


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--threshold", type=int, default=2000,
                        help="Allowable overlap (in nucleotides) (default: %(default)s)")
    parser.add_argument("files", type=str, nargs="*", default=[],
                        help="Specific files to compare (default: all files in 'data' and 'pending')")

    args = parser.parse_args()
    if run_all(args.files, args.threshold):
        sys.exit(1)
