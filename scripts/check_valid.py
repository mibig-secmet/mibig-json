#!/usr/bin/env python3

from collections import defaultdict
import glob
import json
import os
import sys
from typing import Any, Dict, List

from jsonschema import validate, ValidationError

with open("schema.json") as schema_handle:
    schema = json.load(schema_handle)


def check_coordinates(data: Dict[str, Any], prefix: str) -> bool:
    start = data["cluster"]["loci"].get("start_coord", 1)
    end = data["cluster"]["loci"].get("end_coord", start + 10)
    if start >= end:
        print(f"{prefix}invalid coordinates: end ({end}) before start ({start})")
        return False
    return True


def check_gene_naming(data: Dict[str, Any], prefix: str) -> bool:
    invalid = {
        "no_accession",
        "unknown",
    }
    ids = set()
    for gene in data["cluster"].get("genes", {}).get("annotations", {}):
        ids.add(gene["id"])
        name = gene.get("name")
        if name:
            ids.add(name)
    nrp = data["cluster"].get("nrp", {})
    if nrp:
        for nrps_gene in nrp.get("nrps_genes", []):
            ids.add(nrps_gene["gene_id"])
        for thioesterase in nrp.get("thioesterases", []):
            ids.add(thioesterase["gene"])
    for synthase in data["cluster"].get("polyketide", {}).get("synthases", []):
        for gene in synthase.get("genes", []):
            ids.add(gene)
        for module in synthase.get("modules", []):
            for gene in module.get("genes", []):
                ids.add(gene)
    terpene = data["cluster"].get("terpene", {})
    if terpene:
        ids.update(terpene.get("terpene_prenyltransferases", []))
        ids.update(terpene.get("terpene_synth_cycl", []))
    bad_names = {name.lower() for name in ids}.intersection(invalid)
    if bad_names:
        print(f"{prefix}invalid gene identifiers: {', '.join(list(bad_names))}")
        return False
    return True


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


def check_item_duplication(data: Dict[str, Any], prefix: str) -> bool:
    def check_list(items) -> bool:
        if items and isinstance(items[0], str):
            return len(items) == len(set(items))
        return True

    for key, val in data.items():
        if isinstance(val, list) and not check_list(val):
            print(f"{prefix}contains duplicated items in '{key}': {val}")
            return False
        if isinstance(val, dict):
            if not check_item_duplication(val, prefix):
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


def check_chem_act(data: Dict[str, Any], prefix: str) -> bool:
    comp_act: Dict[str, List[str]] = defaultdict(list)
    for compound in data["cluster"].get("compounds", []):
        for activity in compound.get("chem_acts", []):
            if activity["activity"].lower() in ("unknown"):
                comp_act[compound["compound"]].append(activity["activity"])
    if comp_act:
        print(f"{prefix}contains compounds with invalid activities:")
        for name, activities in comp_act.items():
            print(f"\t{name}: {activities}")
        return False
    return True


def check_nrp_opts(data: Dict[str, Any], prefix: str) -> bool:
    if "nrp" in data["cluster"] and "NRP" not in data["cluster"]["biosyn_class"]:
        print(f"{prefix}contains nrp details but doesn't list 'NRP' as biosynthetic class")
        return False
    return True


def check_all() -> bool:
    valid = check_multiple(glob.glob(os.path.join("data", "*.json")), prefix="data/")
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
            check_coordinates,
            check_gene_duplication,
            check_gene_naming,
            check_item_duplication,
            check_kr_stereochem,
            check_pks_module_duplication,
            check_chem_act,
            check_nrp_opts,
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
