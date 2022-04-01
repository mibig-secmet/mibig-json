#!/usr/bin/env python3

import glob
import json
import os
import sys
from typing import List

from jsonschema import validate, ValidationError

with open("schema.json") as schema_handle:
    schema = json.load(schema_handle)


def check_all() -> bool:
    valid = True
    for dir_name in ["data", "pending"]:
        valid = check_multiple(glob.glob(os.path.join(dir_name, "*.json"))) and valid
    return valid


def check_multiple(files: List[str]) -> bool:
    valid = True
    for file in files:
        valid = check_single(file) and valid
    return valid


def check_single(file: str) -> bool:
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
    return True


if __name__ == "__main__":
    if len(sys.argv) < 2:
        if not check_all():
            sys.exit(1)
    else:
        if not check_multiple(sys.argv[1:]):
            sys.exit(1)
