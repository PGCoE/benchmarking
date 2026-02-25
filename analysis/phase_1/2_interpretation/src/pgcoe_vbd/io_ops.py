from pathlib import Path
from typing import List
import pandas as pd
import streamlit as st
import screed

SAMPLESHEET_COLS = [
    "sample", "species", "segment", 
    "source", "assembly"
]

def load_csv(path: Path, required_cols: List = []) -> pd.DataFrame:    

    df = pd.read_csv(path)

    # Validate required columns are present
    missing = set(required_cols) - set(df.columns)
    if missing:
        st.error(f"{path}: missing required columns: {sorted(missing)}")

    if df.empty:
        raise ValueError("No input rows found.")

    return df


def load_one_sequence(path: str) -> str:
    """
    Load ONLY the first contig/record from a FASTA/FASTQ using screed.
    """

    with screed.open(path) as records:
        for rec in records:
            return rec["sequence"]   # first contig only

    raise ValueError(f"No sequences found in: {path}")


