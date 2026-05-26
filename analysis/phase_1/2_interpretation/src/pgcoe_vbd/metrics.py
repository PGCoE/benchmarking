#!/usr/bin/env python3

import subprocess
import sys
from typing import List, Dict, Any, Tuple
import os
import screed
import numpy as np
import streamlit as st
from Bio.Align import PairwiseAligner
from collections import defaultdict
import itertools
import math
import pandas as pd

# IUPAC ambiguity codes for nucleotide matching
IUPAC_CODES = {
    'A': {'A'},
    'C': {'C'},
    'G': {'G'},
    'T': {'T'},
    'U': {'T'},  # Treat U as T
    'R': {'A', 'G'},
    'Y': {'C', 'T'},
    'S': {'G', 'C'},
    'W': {'A', 'T'},
    'K': {'G', 'T'},
    'M': {'A', 'C'},
    'B': {'C', 'G', 'T'},
    'D': {'A', 'G', 'T'},
    'H': {'A', 'C', 'T'},
    'V': {'A', 'C', 'G'}
}

def calculate_sequence_metrics(seq1, seq2) -> Dict[str, Any]:
    """
    Calculate similarity metrics between two aligned sequences using NumPy arrays.
    Does NOT compute match_iupac (IUPAC ambiguities are treated as mismatches for match/sub purposes,
    except invalid characters are still tracked).
    """
    s1 = np.array(list(seq1.sequence.upper()))
    s2 = np.array(list(seq2.sequence.upper()))

    if s1.shape != s2.shape:
        sys.exit("Error: Sequences are not the same length!")

    metrics: Dict[str, Any] = {
        'seq1': seq1.name,
        'seq2': seq2.name,
        'match': 0,
        'ins': 0,
        'del': 0,
        'sub': 0,
        'invalid': 0,
        'valid': 0,
        'valid_s1': 0,
        'invalid_s1': 0,
        'valid_s2': 0,
        'invalid_s2': 0,
        'invalid_char_s1': [],
        'invalid_char_s2': [],
        'del_terminal': 0,
        'del_internal': 0,
        'ins_terminal': 0,
        'ins_internal': 0,
    }

    # Validate sequences against IUPAC codes (still used only for validity)
    valid_bases = np.array(list(IUPAC_CODES.keys()) + ['-'])
    valid_s1 = np.isin(s1, valid_bases)
    valid_s2 = np.isin(s2, valid_bases)

    # Create mask for valid positions (excluding double gaps)
    mask = valid_s1 & valid_s2 & ~((s1 == '-') & (s2 == '-'))

    # Record validation statistics
    metrics['valid_s1'] = int(valid_s1.sum())
    metrics['invalid_s1'] = int((~valid_s1).sum())
    metrics['valid_s2'] = int(valid_s2.sum())
    metrics['invalid_s2'] = int((~valid_s2).sum())
    metrics['valid'] = int(mask.sum())
    metrics['invalid'] = int((~mask).sum())
    metrics['invalid_char_s1'] = list(set(list(s1[~valid_s1])))
    metrics['invalid_char_s2'] = list(set(list(s2[~valid_s2])))

    # Extract valid positions
    b1 = s1[mask]
    b2 = s2[mask]
    length = len(b1)
    position_indices = np.arange(length)

    # Find sequence boundaries (first and last non-gap positions)
    non_gap_b1 = np.flatnonzero(b1 != '-')
    first_b1 = int(non_gap_b1[0]) if non_gap_b1.size > 0 else length
    last_b1  = int(non_gap_b1[-1]) if non_gap_b1.size > 0 else -1

    non_gap_b2 = np.flatnonzero(b2 != '-')
    first_b2 = int(non_gap_b2[0]) if non_gap_b2.size > 0 else length
    last_b2  = int(non_gap_b2[-1]) if non_gap_b2.size > 0 else -1

    # Count exact matches only
    metrics['match'] = int((b1 == b2).sum())

    # Calculate indels
    deletions_mask  = (b2 == '-') & (b1 != '-')
    insertions_mask = (b1 == '-') & (b2 != '-')

    metrics['del'] = int(deletions_mask.sum())
    metrics['ins'] = int(insertions_mask.sum())

    # Terminal vs internal deletions (using seq2 boundaries)
    del_terminal_mask = deletions_mask & ((position_indices < first_b2) | (position_indices > last_b2))
    metrics['del_terminal'] = int(del_terminal_mask.sum())
    metrics['del_internal'] = int(metrics['del'] - metrics['del_terminal'])

    # Terminal vs internal insertions (using seq1 boundaries)
    ins_terminal_mask = insertions_mask & ((position_indices < first_b1) | (position_indices > last_b1))
    metrics['ins_terminal'] = int(ins_terminal_mask.sum())
    metrics['ins_internal'] = int(metrics['ins'] - metrics['ins_terminal'])

    # Substitutions: anything different where neither is a gap
    # (Ambiguous IUPAC codes are treated as normal symbols here; no match_iupac logic.)
    metrics['sub'] = int(((b1 != b2) & (b1 != '-') & (b2 != '-')).sum())

    return metrics


def read_alignment_file(input_file: str) -> List:
    """Read alignment file and return list of sequence records."""
    st.session_state.submsg.info('Reading alignment')
    records = list(screed.open(input_file))
    for record in records:
        if record.name.startswith('_R_'):
            record.name = record.name[3:]
    return records


def run_mafft_alignment(sample_name: str, input_file: str) -> str:
    st.session_state.submsg.info('Aligning sequences with MAFFT')
    output_file = os.path.join(st.session_state.outdir, f'{sample_name}.aligned.fa')

    with open(output_file, 'w') as out:
        subprocess.run(
            ["mafft", "--auto", "--adjustdirection", input_file],
            stdout=out,
            stderr=subprocess.DEVNULL,
            check=True
        )
    return output_file


def create_multi_fasta(sample_name: str, data: List[Dict]) -> Tuple[str, Dict[str, int]]:
    output_file = os.path.join(st.session_state.outdir, f'{sample_name}.multi.fa')
    lengths: Dict[str, int] = {}

    with open(output_file, 'w') as out:
        for row in data:
            source = (row.get('source') or "").strip()
            assembly = (row.get('assembly') or "").strip()
            if not source or not assembly:
                continue

            sequence_count = 0
            for record in screed.open(assembly):
                sequence_count += 1
                if sequence_count > 1:
                    st.error(f'Error: Each file must contain exactly one sequence ({assembly})')
                    return

                out.write(f">{source}\n{record.sequence}\n")
                lengths[source] = len(record.sequence)

    return output_file, lengths


def reverse_complement(seq: str) -> str:
    complement = str.maketrans('ACGTacgt', 'TGCAtgca')
    return seq.translate(complement)[::-1]

def align_pair(seq1, seq2):
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -1

    fwd_alignment = aligner.align(seq1, seq2)
    rev_alignment = aligner.align(reverse_complement(seq1), seq2)

    if rev_alignment[0].score > fwd_alignment[0].score:
        return rev_alignment[0], True  # True = was reoriented
    return fwd_alignment[0], False


def run():
    st.session_state.msg.empty()
    st.session_state.submsg.empty()

    samplesheet = st.session_state.df_samplesheet

    if not st.session_state.outdir.exists():
        st.session_state.outdir.mkdir()
    
    groups = defaultdict(list)
    unique_keys = set()

    for i, record in enumerate(samplesheet.to_dict(orient="records")):
        unique_keys.update(record.keys())

        sample = str(record.get("sample", "")).strip()
        source = str(record.get("source", "")).strip()

        if not sample:
            st.error(f"Record {i} missing 'sample' column")
            return
        elif not source:
            st.error(f"Record {i} missing 'source' column")
            return
        elif source.startswith("_R_"):
            st.error(f"Source cannot start with '_R_'. Please update record {i}: {source}")
            return
        else:
            groups[sample].append(record)

    out: List[Dict[str, Any]] = []

    # Pre-compute total pairs across all samples for global progress
    total_pairs = sum(
        math.comb(len(data), 2) for data in groups.values()
    )
    completed_pairs = 0

    for sample, data in groups.items():

        st.session_state.msg.info(f"Processing {sample}")

        multi_fasta_file, lengths = create_multi_fasta(sample, data)
        aligned_fasta_file = run_mafft_alignment(sample, multi_fasta_file)
        aligned_records = read_alignment_file(aligned_fasta_file)

        data_map = {rec['source']: rec for rec in data}

        for seq1, seq2 in itertools.combinations(aligned_records, 2):

            completed_pairs += 1
            st.session_state.progress.progress(completed_pairs / total_pairs, text=f"Calculating metrics ({completed_pairs} / {total_pairs})")
            st.session_state.submsg.info(f"{seq1.name} vs {seq2.name}")

            pair = calculate_sequence_metrics(seq1, seq2)
            pair['sample'] = sample
            pair['length_s1'] = lengths.get(seq1.name)
            pair['length_s2'] = lengths.get(seq2.name)

            for k in unique_keys:
                pair[k + '_s1'] = data_map[seq1.name].get(k, None)
                pair[k + '_s2'] = data_map[seq2.name].get(k, None)

            out.append(pair)

    st.session_state.df = pd.DataFrame(out)
    st.session_state.msg.empty()
    st.session_state.submsg.empty()
    st.session_state.progress.empty()