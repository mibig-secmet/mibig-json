#!/bin/bash

set -o errexit
set -o pipefail

if [ $# -lt 1 ]; then
    echo "Usage: $0 <version> [mibig_data_dir]"
    exit 2
fi

VERSION="$1"
DATA_DIR="${2:-data}"
OUT_DIR="mibig_json_${VERSION}"

mkdir "$OUT_DIR"

for FILE in "$DATA_DIR"/*.json; do
    STATUS=$(jq .cluster.status "$FILE" | tr -d '"')
    if [ -n "$STATUS" ] && [ "$STATUS" = "active" ]; then
        cp "$FILE" "$OUT_DIR"
    fi
done

tar -czf "${OUT_DIR}.tar.gz" "$OUT_DIR"
