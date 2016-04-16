#!/bin/bash
set -e

dir="tests/data"
mkdir -p "$dir"

function get_dataset {
    local path="$1"
    local filename="$2"
    local tmp_data="$dir/$2"
    local url="https://$path/$filename"
    echo "Downloading $tmp_data from $url..."
    curl -s "$url" > "$tmp_data"
    echo "Uncompressing $tmp_data"
    tar -xzf "$tmp_data" -C "$dir"
}

function get_dataset1 {
    local path="s3.amazonaws.com/s3.aws.paperg.com/wintermute/test-data/2015_12_16_full_staining"
    get_dataset $path "1-1.fcs.tar.gz"
    get_dataset $path "1-1_transformed.txt.tar.gz"
}

function get_dataset2 {
    local path="s3.amazonaws.com/s3.aws.paperg.com/wintermute/test-data/2016_01_09"
    get_dataset $path "2-4.fcs.tar.gz"
    get_dataset $path "2-4_transformed.txt.tar.gz"
}

get_dataset1
get_dataset2
