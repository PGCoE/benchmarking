#!/bin/bash

# samplesheet must use paths that match $PWD

docker run \
    --rm -v "$PWD":"$PWD" \
    -p 8501:8501 \
    public.ecr.aws/o8h2f0o1/pgcoe_vb_phase1_interpretation:1.1