#!/usr/bin/env python3

# gather-pairs.py
# Author: Jared Johnson, jared.johnson@doh.wa.gov

import argparse
import os
import textwrap
import csv
import logging
import sys
from collections import defaultdict

# ----- FUNCTIONS ----- #

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

def exportData(data, filename):
    if len(data) > 0:
        with open(filename, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=data[0].keys())
            writer.writeheader()
            writer.writerows(data)
        logging.info(f'Data saved to {filename}')
    
def main():
    version = "1.0"

    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=str, nargs='+', required=True, help="Samplehseet paths")
    parser.add_argument('--version', action='version', version=version)
    
    args = parser.parse_args()

    start = f"""
    gather-pairs.py v{version}
    """
    print(textwrap.dedent(start), flush=True)

    joined = defaultdict(list)
    for file in args.input:
        with open(file, newline='') as f:
            reader = csv.DictReader(f)
            for row in reader:
                joined[row['sample']].append(row)
    
    singles, pairs = [], []
    for k, v in joined.items():
        if len(v) == 1:
            singles.append(v[0])
        else:
            for i in v:
                pairs.append(i)

    exportData(pairs, f'samplesheet.pairs.csv')
    exportData(singles, f'samplesheet.singles.csv')
            
if __name__ == "__main__":
    main()