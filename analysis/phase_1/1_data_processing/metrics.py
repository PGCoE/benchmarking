#!/usr/bin/env python3

# metrics.py
# Author: Jared Johnson, jared.johnson@doh.wa.gov

import argparse
import csv
import logging
import subprocess
import sys
import textwrap
from collections import defaultdict
from typing import List, Dict, Any, Tuple
import os

import screed
import numpy as np
import itertools

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

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

# --------------------------
# Core Functions
# --------------------------

def load_csv_data(file_path: str, delimiter: str = ',') -> List[Dict[str, str]]:
    """Load CSV data and return list of dictionaries."""
    with open(file_path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter=delimiter))


def calculate_sequence_metrics(seq1, seq2) -> Dict[str, Any]:
    """
    Calculate similarity metrics between two aligned sequences using NumPy arrays.
    
    Args:
        seq1: First aligned sequence record
        seq2: Second aligned sequence record
        
    Returns:
        Dictionary containing alignment metrics
    """
    logging.info(f'Calculating metrics (ref: {seq1.name}, query: {seq2.name})')

    s1 = np.array(list(seq1.sequence.upper()))
    s2 = np.array(list(seq2.sequence.upper()))

    if s1.shape != s2.shape:
        sys.exit("Error: Sequences are not the same length!")

    # Initialize metrics dictionary
    metrics = {
        'match': 0,
        'match_iupac': 0,  # Add to subs when ignoring IUPAC matches
        'ins': 0,
        'del': 0,
        'sub': 0,
        'invalid': 0,
        'valid': 0,
        'valid_s1': 0,
        'invalid_s1': 0,
        'valid_s2': 0,
        'invalid_s2': 0,
        'del_terminal': 0,
        'del_internal': 0,
        'ins_terminal': 0,
        'ins_internal': 0,
    }
    metrics['seq1'] = seq1.name
    metrics['seq2'] = seq2.name

    # Validate sequences against IUPAC codes
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
    last_b1 = int(non_gap_b1[-1]) if non_gap_b1.size > 0 else -1

    non_gap_b2 = np.flatnonzero(b2 != '-')
    first_b2 = int(non_gap_b2[0]) if non_gap_b2.size > 0 else length
    last_b2 = int(non_gap_b2[-1]) if non_gap_b2.size > 0 else -1

    # Count matches
    metrics['match'] = int((b1 == b2).sum())

    # Calculate indels (relative to seq1)
    deletions_mask = (b2 == '-') & (b1 != '-')  # deletion in seq2 vs seq1
    insertions_mask = (b1 == '-') & (b2 != '-')  # insertion in seq2 vs seq1

    metrics['del'] = int(deletions_mask.sum())
    metrics['ins'] = int(insertions_mask.sum())

    # Classify deletions as terminal vs internal (using seq2 boundaries)
    del_terminal_left = deletions_mask & (position_indices < first_b2)
    del_terminal_right = deletions_mask & (position_indices > last_b2)
    del_terminal_count = int(del_terminal_left.sum() + del_terminal_right.sum())
    
    del_terminal_mask = deletions_mask & ((position_indices < first_b2) | (position_indices > last_b2))
    metrics['del_terminal'] = int(del_terminal_mask.sum())
    metrics['del_internal'] = int(deletions_mask.sum() - metrics['del_terminal'])

    # Classify insertions as terminal vs internal (using seq1 boundaries)
    ins_terminal_left = insertions_mask & (position_indices < first_b1)
    ins_terminal_right = insertions_mask & (position_indices > last_b1)
    ins_terminal_count = int(ins_terminal_left.sum() + ins_terminal_right.sum())
    
    ins_terminal_mask = insertions_mask & ((position_indices < first_b1) | (position_indices > last_b1))
    metrics['ins_terminal'] = int(ins_terminal_mask.sum())
    metrics['ins_internal'] = int(insertions_mask.sum() - metrics['ins_terminal'])

    # Calculate substitutions with IUPAC awareness
    difference_mask = (b1 != b2)

    # Handle simple bases (A/T/C/G) separately from IUPAC codes
    single_bases = np.array(list('ATCG-'))
    single_base_mask = (np.isin(b1, single_bases) & 
                       np.isin(b2, single_bases) & 
                       difference_mask)

    # Count simple substitutions (both non-gap)
    simple_substitutions = (single_base_mask & 
                          (b1 != '-') & 
                          (b2 != '-'))
    metrics['sub'] = int(simple_substitutions.sum())

    # Handle IUPAC code differences (both not simple bases and both not gaps)
    iupac_differences = (difference_mask & 
                        ~single_base_mask & 
                        (b1 != '-') & 
                        (b2 != '-'))

    # Evaluate IUPAC matches vs true substitutions
    if np.any(iupac_differences):
        for base1, base2 in zip(b1[iupac_differences], b2[iupac_differences]):
            if IUPAC_CODES[base1] & IUPAC_CODES[base2]:
                metrics['match_iupac'] += 1
            else:
                metrics['sub'] += 1

    return metrics


def read_alignment_file(input_file: str) -> List:
    """Read alignment file and return list of sequence records."""
    logging.info('Reading alignment')
    return list(screed.open(input_file))


def run_mafft_alignment(sample_name: str, input_file: str, outdir: str) -> str:
    """
    Run MAFFT alignment on input sequences.
    
    Args:
        sample_name: Name for output file prefix
        input_file: Path to input FASTA file
        
    Returns:
        Path to aligned FASTA file
    """
    logging.info('Aligning sequences with MAFFT')
    output_file = os.path.join(outdir, f'{sample_name}.aligned.fa')
    
    with open(output_file, 'w') as out:
        subprocess.run(
            ["mafft", "--auto", input_file],
            stdout=out,
            stderr=subprocess.DEVNULL,
            check=True
        )
    return output_file


def create_multi_fasta(sample_name: str, sample_data: List[Dict[str, str]], outdir: str) -> str:
    """
    Combine sequences from multiple files into a multi-FASTA file.
    
    Args:
        sample_name: Name for output file prefix
        sample_data: List of dictionaries containing submitter and assembly info
        
    Returns:
        Path to created multi-FASTA file
    """
    logging.info('Combining sequences into a multi-FASTA')
    output_file = os.path.join(outdir, f'{sample_name}.multi.fa')
    lengths = {}
    with open(output_file, 'w') as out:
        for row in sample_data:
            submitter = row.get('submitter')
            assembly = row.get('assembly')
            
            if not submitter or not assembly:
                continue
                
            sequence_count = 0
            for record in screed.open(assembly):
                sequence_count += 1
                if sequence_count > 1:
                    sys.exit(f'Error: Each file must contain exactly one sequence ({assembly})')
                    
                out.write(f">{submitter}\n{record.sequence}\n")
                lengths[submitter] = len(record.sequence)
                
    return output_file, lengths


def save_metrics_to_csv(metrics_data: List[Dict[str, Any]], filename: str = "metrics.csv") -> None:
    """
    Save list of metrics to CSV file.
    
    Args:
        metrics_data: List of metric dictionaries
        filename: Output filename
    """
    if not metrics_data:
        logging.warning(f"No data to save for {filename}")
        return
        
    fieldnames = list(metrics_data[0].keys())
    
    with open(filename, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(metrics_data)
        
    logging.info(f'Data saved to {filename}')


# --------------------------
# Main
# --------------------------

def main():
    """Main execution function."""
    version = "1.1"
    parser = argparse.ArgumentParser(
        description="Calculate sequence similarity metrics between aligned sequences."
    )
    parser.add_argument("--input", type=str, help="Path to samplesheet file")
    parser.add_argument("--outdir", type=str, default='./', help="Output directory")
    parser.add_argument("--test", action="store_true", 
                       help="Process only first sample then stop.")
    parser.add_argument('--version', action='version', 
                       version=f'%(prog)s {version}')
    args = parser.parse_args()

    print(textwrap.dedent(f"""\
        metrics.py v{version}
        ----------------------
    """), flush=True)
    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)

    # Load samplesheet and group by sample
    samplesheet = load_csv_data(args.input)
    sample_groups = defaultdict(list)
    
    for record in samplesheet:
        sample = record.get('sample')
        if sample:
            sample_groups[sample].append(record)

    # Process each sample group
    for count, (sample_name, sample_data) in enumerate(sample_groups.items(), 1):
        print(f"{sample_name} has {len(sample_data)} assemblies")
        
        # Create multi-FASTA and align sequences
        multi_fasta_file, lengths = create_multi_fasta(sample_name, sample_data, outdir)
        aligned_fasta_file = run_mafft_alignment(sample_name, multi_fasta_file, outdir)
        aligned_records = read_alignment_file(aligned_fasta_file)

        # Calculate all pairwise metrics
        pairwise_metrics = []
        for seq1, seq2 in itertools.combinations(aligned_records, 2):
            metrics = calculate_sequence_metrics(seq1, seq2)
            metrics['sample'] = sample_name
            metrics['length_s1'] = lengths[seq1.name]
            metrics['length_s2'] = lengths[seq2.name]
            pairwise_metrics.append(metrics)

        # Save results
        if pairwise_metrics:
            metrics_file = os.path.join(outdir, f"{sample_name}_metrics.csv")
            save_metrics_to_csv(pairwise_metrics, metrics_file)
        else:
            logging.warning(f"No pairwise comparisons computed for {sample_name}")
            
        # Exit early if in test mode
        if args.test and count == 1:
            sys.exit(0)


if __name__ == "__main__":
    main()