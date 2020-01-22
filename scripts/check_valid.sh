#!/bin/bash

set -o errexit
set -o pipefail

FAILED=()

for changed in data/*.json pending/*.json; do
    jsonschema -i $changed schema.json || FAILED+=($changed)
done

if [ ${#FAILED} -gt 0 ]; then
    echo "Validation failed on the following files: ${FAILED[@]}"
    exit 1
fi
