#!/usr/bin/env python3
"""Update the status of a file"""


import argparse
import json


def run(args: argparse.Namespace) -> None:
    for json_file in args.file:
        with open(json_file, 'r', encoding="utf-8") as handle:
            entry = json.load(handle)

        cluster = entry['cluster']

        old_status = cluster['status'] if 'status' in cluster else "unknown"

        cluster['status'] = args.status
        if args.status == "retired":
            cluster['retirement_reasons'] = [args.reason]
            if args.see_also:
                cluster['see_also'] = args.see_also
        else:
            if old_status == "retired":
                del cluster['retirement_reasons']
                if 'see_also' in cluster:
                    del cluster['see_also']

        if args.comment:
            entry['comments'] = args.comment

        with open(json_file, 'w', encoding="utf-8") as handle:
            json.dump(entry, handle, indent=4, separators=(
                ',', ': '), sort_keys=True, ensure_ascii=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("file", nargs="+",
                        help="File(s) to update status in")
    parser.add_argument("--status", "-s", default="active",
                        choices=("pending", "active", "retired"),
                        help="Status to set for the entry (default: %(default)s).")
    parser.add_argument("--reason", "-r",
                        help="Reason for retiring a file (ignored for other status values)")
    parser.add_argument("--see-also", "-S", action="append",
                        help="Non-retired records to reference (ignored for other status values)")
    parser.add_argument("--comment", "-c",
                        help="Comment to add to entry")
    args = parser.parse_args()
    if args.status == "retired":
        if not args.reason:
            parser.error("Retiring an entry needs a --reason.")

    run(args)


if __name__ == "__main__":
    main()
