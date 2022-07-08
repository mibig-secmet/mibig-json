#!/usr/bin/env python3
"""Try and make sure compound names are following the same rules."""


import argparse
import json
import os
from typing import List


DATA_DIR = os.path.join(os.path.abspath(os.path.dirname(os.path.dirname(__file__))), "data")


def run(args: argparse.Namespace) -> None:
    json_files: List[str] = []
    for name in os.listdir(args.data_dir):
        if name.startswith("BGC") and name.endswith(".json"):
            json_files.append(os.path.join(args.data_dir, name))

    for json_file in json_files:
        rename_compounds(json_file)


def rename_compounds(json_file: str) -> None:
    with open(json_file, 'r', encoding="utf-8") as handle:
        entry = json.load(handle)

    for compound_data in entry['cluster']['compounds']:
        compound_data['compound'] = fix_name(compound_data['compound'])

    with open(json_file, 'w', encoding="utf-8") as handle:
        json.dump(entry, handle, indent=4, separators=(',', ': '), sort_keys=True,
                  ensure_ascii=False)


def fix_name(name: str) -> str:
    parts = name.strip().split(" ")
    parts[0] = fix_string(parts[0])
    for i, part in enumerate(parts):
        parts[i] = fix_greek_letters(part)
    new_name = " ".join(parts)
    return new_name


def fix_greek_letters(name: str) -> str:
    parts = name.split('-')
    for i, part in enumerate(parts):
        part = replace_greek(part)
        if "," in part:
            comma_parts = part.split(',')
            for j, comma_part in enumerate(comma_parts):
                comma_parts[j] = replace_greek(comma_part)
            part = ",".join(comma_parts)
        parts[i] = part
    return "-".join(parts)


def replace_greek(symbol: str) -> str:
    if symbol == "alpha":
        return "α"
    elif symbol == "beta":
        return "β"
    elif symbol == "gamma":
        return "γ"
    elif symbol == "delta":
        return "δ"
    elif symbol == "epsilon":
        return "ε"
    elif symbol == "zeta":
        return "ζ"
    elif symbol == "eta":
        return "η"
    elif symbol == "theta":
        return "θ"
    elif symbol == "iota":
        return "ι"
    elif symbol == "kappa":
        return "κ"
    elif symbol == "lambda":
        return "λ"
    elif symbol == "mu":
        return "μ"
    elif symbol == "nu":
        return "ν"
    elif symbol == "xi":
        return "ξ"
    elif symbol == "omicron":
        return "ο"
    elif symbol == "pi":
        return "π"
    elif symbol == "rho":
        return "ρ"
    elif symbol == "sigma":
        return "σ"
    elif symbol == "tau":
        return "τ"
    elif symbol == "upsilon":
        return "υ"
    elif symbol == "phi":
        return "φ"
    elif symbol == "chi":
        return "χ"
    elif symbol == "psi":
        return "ψ"
    elif symbol == "omega":
        return "ω"
    return symbol


SKIPPABLE_PREFIX_CHARS = set("0123456789,")


def fix_string(name: str) -> str:
    if name.upper() == name:
        return name
    if len(name) <= 3:
        if name.lower() == "iso":
            return "iso"
        return name
    if 3 <= len(name) <= 5 and name[0].isupper() and name[1].islower() and name[-1].isupper():
        return name
    if 3 <= len(name) <= 5 and name[0].isupper() and name[1].islower() and name[-1].isdigit() \
            and name[-2].isupper():
        return name
    if name.startswith("CDA"):
        return name
    if "-" in name:
        parts = name.split("-")
        for i, part in enumerate(parts):
            if not set(part).difference(SKIPPABLE_PREFIX_CHARS):
                continue
            if part in ('L', 'D'):
                continue
            parts[i] = fix_string(part)
            break

        return "-".join(parts)
    return name.lower()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--data-dir', default=DATA_DIR,
                        help="Directory containing the MIBiG JSON files (default: %(default)s).")

    args = parser.parse_args()

    run(args)


if __name__ == "__main__":
    main()
