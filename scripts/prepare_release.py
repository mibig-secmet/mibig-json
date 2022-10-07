#!/usr/bin/env python3

from argparse import ArgumentParser
import glob
import json
import os
from typing import List

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def main() -> None:
    parser = ArgumentParser()
    parser.add_argument("version",
                        help="Version number to replace 'next' version entries by.")
    parser.add_argument("files", nargs="*",
                        help="File(s) to update the version number on")

    args = parser.parse_args()

    if not args.files:
        args.files = sorted(glob.glob(os.path.join(BASE_DIR, "data", "*.json")))


def run(files: List[str], version: str) -> None:
    version_parts = version.split(".")
    for part in version_parts:
        assert part.isnumeric(), f"invalid version specifier {version}"

    for json_file in files:
        with open(json_file, "r", encoding="utf-8") as handle:
            data = json.load(handle)

        last_log = data["changelog"][-1]
        if last_log["version"] == "next":
            last_log["version"] = version

        for entry in data["changelog"]:
            assert entry["version"] != "next", f"Invalid changelog version in {json_file}"

        with open(json_file, "w", encoding="utf-8") as handle:
            json.dump(data, handle, indent=4, separators=(',', ': '),
                      sort_keys=True, ensure_ascii=False)


if __name__ == "__main__":
    main()
