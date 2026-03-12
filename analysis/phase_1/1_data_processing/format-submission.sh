#!/usr/bin/env bash
set -e
shopt -s nullglob  # if a glob has no matches, it expands to nothing (skip cleanly)

mkdir -p prepared_entries

for d in data_formatted/*; do
  [ -d "$d" ] || continue
  echo -e "\nStarting $d"

  # pick the first .xlsx in the dataset (assumes exactly one per folder)
  meta=( "$d"/*.xlsx )
  [ ${#meta[@]} -gt 0 ] || { echo "  No .xlsx in $d, skipping"; continue; }

  docker run --rm \
    -v "$PWD":/"$PWD" \
    public.ecr.aws/o8h2f0o1/pgcoe_vb_phase1_data_processing:1.1 \
    ./prepare-entry.py \
      --meta "$PWD/${meta[0]}" \
      --fasta "$PWD/$d"/assemblies/* \
      --outdir "$PWD/prepared_entries/" \
      --workflow-map "$PWD"/workflow_map.csv
done

docker run --rm \
    -v "$PWD":"$PWD" \
    public.ecr.aws/o8h2f0o1/pgcoe_vb_phase1_data_processing:1.1 \
    ./gather-pairs.py \
      --input "$PWD"/prepared_entries/*.csv \
      --outdir "$PWD"