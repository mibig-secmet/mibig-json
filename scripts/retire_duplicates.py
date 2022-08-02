#!/usr/bin/env python3

from argparse import ArgumentParser
import glob
import json
import os
import re
import subprocess
from typing import List


RETIRED = re.compile(r"^Retired as.*duplicate of (BGC.*)$")


def main() -> None:
    parser = ArgumentParser()
    parser.add_argument("files", nargs="*",
                        help="Files to check and potentially retire for being duplicates")
    args = parser.parse_args()

    if not args.files:
        args.files = sorted(glob.glob(os.path.join("retired", "*.json")))

    run(args.files)


def run(files: List[str]) -> None:
    for file in files:
        with open(file, 'r', encoding="utf-8") as handle:
            data = json.load(handle)

        found = False
        dup = ""
        last_log = data["changelog"][-1]
        for comment in last_log["comments"]:
            match = RETIRED.search(comment)
            if not match:
                continue
            dup = match.group(1)
            found = True
            break
        if not found:
            continue

        cmd = ["retire", file, dup]

        subprocess.run(cmd)


if __name__ == "__main__":
    main()
