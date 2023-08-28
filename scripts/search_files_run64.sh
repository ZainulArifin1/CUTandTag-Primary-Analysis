#!/bin/bash
file_list=$(find "$(realpath ../data/run64/equalcellno/)" -type f -name '*.fastq.gz')

for file in $file_list; do
    base_name=$(basename "$file")
    echo "$base_name: $file"
done
