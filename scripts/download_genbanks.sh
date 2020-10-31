#!/bin/bash

# Download GenBank files for MIBiG entries
# Requires:
# * jq
# * ncbi-acc-download
# * sed


set -o errexit
set -o pipefail
set -o errexit

DEFAULT_JSONS_DIR="$(dirname $(dirname $0))/data"
DEFAULT_GENBANKS_DIR="$(dirname $(dirname $0))/genbanks"
readonly VERSION="1.0.0"

usage() {
    echo "Download the GenBank files related to the MIBiG entries"
    echo "Usage: $0 [options]"
    echo "Options:"
    echo "    -a <key>  | --api <key>          NCBI Entrez API key to use for downloading"
    echo "    -j <path> | --jsons <path>       Path to directory containing MIBiG JSON files (default: ${DEFAULT_JSONS_DIR})"
    echo "    -g <path> | --genbanks <path>    Path to directory to hold the downloaded GenBank files (default: ${DEFAULT_GENBANKS_DIR})"
    echo "    -V | --version                   Print the script's version"
    echo "    -h | --help                      Print this help text"
}

download() {
    ACC=$1
    GENBANKS_DIR=$2
    EXTRA_PARAMS=""
    if [ ! -z ${API_KEY} ]; then
        EXTRA_PARAMS="--api-key ${API_KEY}"
    fi
    pushd ${GENBANKS_DIR} > /dev/null
    if [ ! -e "${ACC}.gbk" ]; then
        ncbi-acc-download --recursive $ACC $EXTRA_PARAMS
    fi
    popd > /dev/null
}

getAcc() {
    JSON_FILE=$1
    jq '.cluster.loci.accession' $JSON_FILE | sed -e 's/"//g'
}

run() {
    JSONS_DIR=$1
    GENBANKS_DIR=$2
    for JSON_FILE in ${JSONS_DIR}/*.json; do
        ACC=$(getAcc ${JSON_FILE})
        BGC_ID=$(basename ${JSON_FILE/\.json/})
        echo "$BGC_ID: $ACC" 1>&2
        download $ACC ${GENBANKS_DIR}
    done
}

while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -h|--help)
            usage
            exit 0
            ;;
        -j|--jsons)
            JSONS_DIR="$2"
            shift; shift
            ;;
        -g|--genbanks)
            GENBANKS_DIR="$2"
            shift; shift
            ;;
        -V|--version)
            echo $VERSION
            exit 0
            ;;
        -a|--api)
            API_KEY="$2"
            shift; shift
            ;;
        *)
            usage
            exit 2
            ;;
    esac
done

: ${JSONS_DIR:=${DEFAULT_JSONS_DIR}}
: ${GENBANKS_DIR:=${DEFAULT_GENBANKS_DIR}}

if [[ ! -d ${GENBANKS_DIR} ]]; then
    mkdir -p ${GENBANKS_DIR}
fi

run ${JSONS_DIR} ${GENBANKS_DIR}
