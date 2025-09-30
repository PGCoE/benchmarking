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
    -v "$PWD":/data \
    -w /data \
    public.ecr.aws/o8h2f0o1/pgcoe_benchmark_phase1:1.0 \
    ./prepare-entry.py \
      --meta "${meta[0]}" \
      --fasta "$d"/assemblies/* \
      --outdir /data/prepared_entries/
done

docker run --rm \
    -v "$PWD":/data \
    -w /data \
    public.ecr.aws/o8h2f0o1/pgcoe_benchmark_phase1:1.0 \
    ./gather-pairs.py \
      --input prepared_entries/*.csv