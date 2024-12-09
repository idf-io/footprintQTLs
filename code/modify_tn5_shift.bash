#!/usr/bin/env bash

set -euo pipefail

DATA_DIR="$1"
OUT_DIR="$2"
plus_shift_og="$3"
minus_shift_og="$4"

plus_shift=$(( plus_shift_og - 4))
minus_shift=$(( minus_shift_og + 5))

find "$DATA_DIR" -maxdepth 1 -mindepth 1 -type f -name "*sorted.tsv.gz" | while IFS= read -r f; do
	
bgzip -dc "$f" | grep -v "^#" | awk -v ps="$plus_shift" -v ms="$minus_shift" -F'\t' 'BEGIN { OFS="\t" } { $2=$2+ps; $3=$3+ms; print }' - > "${OUT_DIR}/$(basename "${f%_sorted.tsv.gz}")_sorted_ps${plus_shift_og}_ms${minus_shift_og#-}.tsv"

done
