#!/usr/bin/env python3

import argparse
from datetime import datetime
import json
import os
from typing import Any, Dict, List


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--release", default="next",
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
    new = False

    if previous["version"] != release:
        entry = {
            "comments": [],
            "contributors": [],
            "updated_at": [],
            "version": release,
        }
        new = True
    else:
        entry = previous

    assert len(comments) == len(contributors)
    for comment, contributor in zip(comments, contributors):
        entry["comments"].append(comment)
        entry["contributors"].append(contributor)
        entry["updated_at"].append(datetime.now().strftime("%Y-%m-%dT%H:%M:%S%:z"))

    if new:
        record["changelog"].append(entry)

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
