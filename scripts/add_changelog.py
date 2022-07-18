#!/usr/bin/env python3

import argparse
import json
import os
from typing import Any, Dict, List


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--release", default="3.0",
                        help="Set the release version of the changelog entry to add (default: %(default)s).")
    parser.add_argument("-c", "--comment", action="append",
                        help="Comment to add to changelog entry (can be used multiple times)")
    parser.add_argument("-C", "--contributor", action="append",
                        help="Contributor ID of a contributor to this change (can be used multiple times)")
    parser.add_argument("-a", "--amend", action="store_true", default=False,
                        help="Fix the last changelog entry instead of adding a new entry")
    parser.add_argument("json_file", help="Filename of json file to add changelog entry to.")
    args = parser.parse_args()
    if not args.contributor and os.getenv("MIBIG_ID"):
        args.contributor = [os.getenv("MIBIG_ID")]
    run(args)


def run(args: argparse.Namespace) -> None:
    with open(args.json_file, 'r', encoding="utf-8") as handle:
        record = json.load(handle)

    comments = [comment.strip() for comment in args.comment]
    if args.amend:
        amend_changelog(record, args.release, comments, args.contributor)
    else:
        add_changelog(record, args.release, comments, args.contributor)

    with open(args.json_file, 'w', encoding="utf-8") as handle:
        json.dump(record, handle, indent=4, separators=(',', ': '), sort_keys=True, ensure_ascii=False)


def add_changelog(record: Dict[str, Any], release: str, comments: List[str], contributors: List[str]) -> None:

    if not contributors:
        contributors = ["AAAAAAAAAAAAAAAAAAAAAAAA"]

    previous = record["changelog"][-1]

    if previous["version"] != release:
        new_entry = {
            "comments": comments,
            "contributors": contributors,
            "version": release,
        }
        record["changelog"].append(new_entry)
        return

    for comment in comments:
        if comment not in previous["comments"]:
            previous["comments"].append(comment)
    for contributor in contributors:
        if contributor not in previous["contributors"]:
            previous["contributors"].append(contributor)

def amend_changelog(record: Dict[str, Any], release: str, comments: List[str], contributors: List[str]) -> None:
    entry = record["changelog"][-1]
    entry["version"] = release
    if comments:
        if "comments" not in entry or not entry["comments"]:
            entry["comments"] = comments
        else:
            entry["comments"].extend(comments)

    if contributors:
        if "contributors" not in entry or not entry["contributors"]:
            entry["contributors"] = contributors
        else:
            entry["contributors"].extend(contributors)


if __name__ == "__main__":
    main()
