#!/usr/bin/env python3
"""Try and make sure compound names are following the same rules."""


import argparse
import json
import os
from typing import List


DATA_DIR = os.path.join(os.path.abspath(os.path.dirname(os.path.dirname(__file__))), "data")
VERSION = "2.1"
COMMENTS = ["Remove duplicated citations"]
CONTRIBUTORS = ["AAAAAAAAAAAAAAAAAAAAAAAA"]


def run(args: argparse.Namespace) -> None:
    json_files: List[str] = []
    for name in os.listdir(args.data_dir):
        if name.startswith("BGC") and name.endswith(".json"):
            json_files.append(os.path.join(args.data_dir, name))

    for json_file in json_files:
        fix_publcations(json_file, args)


def fix_publcations(json_file: str, args: argparse.Namespace) -> None:
    with open(json_file, 'r', encoding="utf-8") as handle:
        entry = json.load(handle)

    changed = False

    publications = entry['cluster']['publications']
    old_len = len(publications)
    unique_pubs = {}
    for pub in publications:
        unique_pubs[pub] = None

    if len(unique_pubs) != old_len:
        changed = True
    entry['cluster']['publications'] = list(unique_pubs.keys())

    if changed:
        log = entry["changelog"][-1]
        if log['version'] == args.release:
            if "comments" not in log or not log["comments"]:
                log["comments"] = COMMENTS
            else:
                log["comments"].extend(COMMENTS)

            if "contributors" not in log or not log["contributors"]:
                log["contributors"] = args.contributor
            else:
                log["contributors"].extend(args.contributor)
            pass
        else:
            log = {
                "comments": COMMENTS,
                "contributors": args.contributor,
                "version": args.release,
            }

            entry["changelog"].append(log)


    with open(json_file, 'w', encoding="utf-8") as handle:
        json.dump(entry, handle, indent=4, separators=(',', ': '), sort_keys=True, ensure_ascii=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--data-dir', default=DATA_DIR,
                        help="Directory containing the MIBiG JSON files (default: %(default)s).")
    parser.add_argument('-C',  '--contributor', action="append",
                        help="Contributor ID of a contributor to this change (can be used multiple times)")
    parser.add_argument("-r", "--release", default=VERSION,
                        help="Set the release version of the changelog entry to add (default: %(default)s).")
    args = parser.parse_args()

    if not args.contributor:
        args.contributor = CONTRIBUTORS

    run(args)

if __name__ == "__main__":
    main()
