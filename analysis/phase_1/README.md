# Overview
This document describes the workflow for analyzing the **Phase 1 viral benchmarking data**.

# Protocol

> **Note:** Step 1 requires substantial manual curation. If you wish to skip this, a pre-formatted dataset is available on the benchmarking [submission portal](https://benchmarking.nepgcoe.org) (`data_formatted.tar.gz`).

---

## Step 1. Download the raw data
1. Navigate to the [submission portal](https://benchmarking.nepgcoe.org) and log in with your credentials.  
2. Download each lab’s submission directory.  
3. Confirm that each directory contains a properly formatted `Data Capture Table.xlsx` with all **required** fields completed.

---

## Step 2. Prepare the data for processing
> The `prepare-entry.py` script can help identify malformed data.

1. Verify that each `Data Capture Table.xlsx` is complete.  
2. Organize submissions into lab-specific directories:

   ```text
   data_formatted/
   └── lab_1/
       ├── Data Capture Table.xlsx
       └── assemblies/
           ├── sample01_ref01_lab1.fa
           └── ...
   ```

   > Only one contig is allowed per FASTA file.

3. Run the formatting script:

   ```bash
   bash format-submissions.sh
   ```

   This will:  
   - Save individual samplesheets to `prepared_entries/` using `prepare-entry.py`.  
   - Combine samples with ≥ 2 submissions into `samplesheet.pairs.csv` via `gather-pairs.py`.  
   - Write singletons to `samplesheet.singles.csv`.

---

## Step 3. Calculate pairwise metrics
Run the benchmarking metrics with:

```bash
docker run --rm -v "$PWD":/data public.ecr.aws/o8h2f0o1/pgcoe_benchmark_phase1:1.0 \
    ./metrics.py \
    --input /data/samplesheet.pairs.csv \
    --outdir /data/metrics
```

> This step typically completes in just a few minutes.

---

## Step 4. Visualize the results
Launch the Streamlit dashboard:

```bash
docker run --rm -p 8501:8501 -v "$PWD":/data public.ecr.aws/o8h2f0o1/pgcoe_benchmark_phase1:1.0 \
    streamlit run /app/dashboard.py \
        --server.address=0.0.0.0 \
        --server.port=8501 \
        --server.headless=true
```

Then open [http://localhost:8501](http://localhost:8501) in your browser to explore the benchmarking results.
