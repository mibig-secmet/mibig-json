#!/usr/bin/env python3

from argparse import ArgumentParser
import glob
import json
import os
from typing import Any, Dict, List


def main() -> None:
    parser = ArgumentParser()
    parser.add_argument("files", nargs="*",
                        help="File(s) to fix")
    parser.add_argument("--contributors", "-c", action="append",
                        help="Contributor to add (can be specified multiple times)")
    args = parser.parse_args()

    if not args.files:
        args.files = sorted(glob.glob(os.path.join("retired", "*.json")))

    if not args.contributors:
        args.contributors = ["AAAAAAAAAAAAAAAAAAAAAAAA"]

    run(args.files, args.contributors)


def run(files: List[str], contributors: List[str]) -> None:
    for file in files:
        with open(file, 'r', encoding="utf-8") as handle:
            data = json.load(handle)

        fix_contributors(contributors, data)
        remove_empty_lists(data)

        with open(file, 'w', encoding="utf-8") as handle:
            json.dump(data, handle, indent=4, separators=(',', ': '),
                      sort_keys=True, ensure_ascii=False)


def fix_contributors(contributors: List[str], data: Dict[str, Any]) -> None:
    for entry in data["changelog"]:
        if 'contributors' not in entry:
            entry['contributors'] = contributors


def remove_empty_lists(data: Dict[str, Any]) -> None:
    keys_to_remove = set()
    for key, item in data.items():
        if isinstance(item, dict):
            if not len(item):
                keys_to_remove.add(key)
            else:
                remove_empty_lists(item)
        elif isinstance(item, list):
            if not len(item):
                keys_to_remove.add(key)
            else:
                for subitem in item:
                    if isinstance(subitem, dict):
                        remove_empty_lists(subitem)
                    else:
                        break
    for key in keys_to_remove:
        del data[key]


if __name__ == "__main__":
    main()
