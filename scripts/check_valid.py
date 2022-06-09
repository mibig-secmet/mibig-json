#!/usr/bin/env python3

import glob
import json
import os
import sys
from typing import Any, Dict, List

from jsonschema import validate, ValidationError

with open("schema.json") as schema_handle:
    schema = json.load(schema_handle)


def check_gene_duplication(data: Dict[str, Any], prefix: str) -> bool:
    ids = []
    names = []
    for gene in data["cluster"].get("genes", {}).get("annotations", {}):
        ids.append(gene["id"])
        name = gene.get("name")
        if name:
            names.append(name)
    id_mismatch = len(set(ids)) != len(ids)
    name_mismatch = len(set(names)) != len(names)
    if id_mismatch or name_mismatch:
        message_chunks = []
        if id_mismatch:
            message_chunks.append("ids")
        if id_mismatch and name_mismatch:
            message_chunks.append("and")
        if name_mismatch:
            message_chunks.append("names")
        print(f"{prefix}{' '.join(message_chunks)} duplicated in gene annotations")
        return False
    return True


def check_kr_stereochem(data: Dict[str, Any], prefix: str) -> bool:
    for synthase in data["cluster"].get("polyketide", {}).get("synthases", []):
        for module in synthase.get("modules", []):
            if "kr_stereochem" in module and "Ketoreductase" not in module["domains"]:
                print(f"{prefix}contains KR stereochemistry for modules without KR domains")
                return False
    return True


def check_pks_module_duplication(data: Dict[str, Any], prefix: str) -> bool:
    for synthase in data["cluster"].get("polyketide", {}).get("synthases", []):
        module_numbers = set()
        duplicates = set()
        for module in synthase.get("modules", []):
            num = module.get("module_number")
            if num is None:
                continue
            if num in module_numbers:
                duplicates.add(num)
            module_numbers.add(num)
        if duplicates:
            print(f"{prefix}contains duplicate PKS module numbers: {sorted(duplicates)}")
            return False
    return True


def check_all() -> bool:
    valid = True
    for dir_name in ["data", "pending"]:
        valid = check_multiple(glob.glob(os.path.join(dir_name, "*.json")), prefix=f"{dir_name}/") and valid
    return valid


def check_multiple(files: List[str], prefix: str = "") -> bool:
    valid = True
    for file in sorted(files):
        valid = check_single(file, prefix=prefix) and valid
    return valid


def check_single(file: str, prefix: str = "") -> bool:
    with open(file) as handle:
        try:
            data = json.load(handle)
        except json.JSONDecodeError as err:
            print(f"{os.path.split(file)[-1]}: {err}")
            return False
    try:
        validate(instance=data, schema=schema)
    except ValidationError as err:
        lines = str(err).splitlines()
        print(f"{os.path.split(file)[-1]}: {lines[-2]} {lines[0]}")
        return False

    valid = True

    if prefix:
        prefix = f"{prefix}{file.rsplit(prefix, 1)[-1]}: "
    try:
        for func in [
            check_gene_duplication,
            check_kr_stereochem,
            check_pks_module_duplication,
        ]:
            valid = func(data, prefix) and valid
    except KeyError as err:
        raise ValueError(f"missing expected field in {prefix}: {err}") from err
    return valid


if __name__ == "__main__":
    if len(sys.argv) < 2:
        try:
            if not check_all():
                sys.exit(1)
        except KeyboardInterrupt:  # hides long JSON decoding stack traces
            sys.exit(1)
    else:
        if not check_multiple(sys.argv[1:]):
            sys.exit(1)
