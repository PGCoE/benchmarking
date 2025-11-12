#!/usr/bin/env python3

"""
prepare-entry.py
Author: Jared Johnson, jared.johnson@doh.wa.gov

Script to process genomic assembly metadata and FASTA files for submission preparation.
Reads Excel metadata, validates FASTA files, and exports structured data.
"""

import argparse
import csv
import gzip
import json
import logging
import os
import re
import shutil
import sys
import textwrap
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

from openpyxl import load_workbook

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

# Constants
VERSION = "1.1"
REQUIRED_COLUMNS = ['SRA ID', 'Assembly Species / Segment', 'Assembly Filename', 'Assembly Generated?', 'Sample ID']
ASSEMBLY_GENERATED_VALUES = {'true', 'yes', '1', 'y'}

def sanitize_string(name: str, replacement: str = "_") -> str:
    if not name:
        return ""
    sanitized = re.sub(r'[<>:"/\\|?*\s]+', replacement, str(name))
    sanitized = sanitized.strip(replacement)
    return sanitized[:255]

def ensure_outdir(outdir: Path) -> None:
    try:
        outdir.mkdir(parents=True, exist_ok=True)
    except Exception as e:
        sys.exit(f"Error: Could not create outdir {outdir}: {e}")

def load_workbook_data(file_path: str) -> Tuple[str, Dict]:
    try:
        wb = load_workbook(file_path, data_only=True)
    except Exception as e:
        logger.error(f"Failed to load workbook: {e}")
        sys.exit(f"Error: Could not load workbook {file_path}")

    if 'Sheet1' not in wb.sheetnames:
        sys.exit(f'Error: "Sheet1" not found. Available sheets: {", ".join(wb.sheetnames)}')

    ws = wb['Sheet1']

    data = []
    for row in ws.iter_rows(values_only=True):
        data.append([i.strip() if isinstance(i, str) else i for i in row])

    if len(data) < 6:
        sys.exit("Error: Workbook does not have enough rows")

    submitter_row_valid = data[3] and len(data[3]) > 1 and 'Submitter Name' in str(data[3][1] or '')
    sample_id_row_valid = data[5] and len(data[5]) > 1 and 'Sample ID' in str(data[5][1] or '')

    if not submitter_row_valid or not sample_id_row_valid:
        sys.exit('Error: Invalid workbook structure. Expected "Submitter Name" in row 4, "Sample ID" in row 6')

    submitter = data[4][1] if data[4] and len(data[4]) > 1 else ""
    column_names = data[5] if data[5] else []

    if not submitter:
        sys.exit("Error: Submitter name not found")

    try:
        column_indices = {col: column_names.index(col) for col in REQUIRED_COLUMNS}
    except ValueError:
        missing_cols = [col for col in REQUIRED_COLUMNS if col not in column_names]
        sys.exit(f"Error: Missing required columns: {', '.join(missing_cols)}")

    parsed_data = defaultdict(dict)

    for row_idx, row in enumerate(data[6:], start=7):
        if not row or len(row) <= max(column_indices.values()):
            continue
        try:
            sra_id = row[column_indices['SRA ID']]
            segment = row[column_indices['Assembly Species / Segment']]
            sample_id = row[column_indices['Sample ID']]
            assembly_file = row[column_indices['Assembly Filename']]
            assembly_generated = str(row[column_indices['Assembly Generated?']] or '').lower().strip()

            if not sra_id or not segment:
                continue

            if assembly_generated in ASSEMBLY_GENERATED_VALUES:
                parsed_data[str(sra_id)][str(segment)] = {
                    'name': sample_id,
                    'file': assembly_file,   # may be basename or path; resolved later
                    'submitter': submitter
                }

        except (IndexError, KeyError) as e:
            logger.warning(f"Skipping row {row_idx} due to parsing error: {e}")
            continue

    if not parsed_data:
        sys.exit("Error: No valid assembly data found in workbook")

    return submitter, dict(parsed_data)

def stage_gzip_to_outdir(src: Path, dst_dir: Path) -> Path:
    """
    Ensure a gzipped copy lives in dst_dir. If src is already .gz, copy.
    Otherwise gzip to <dst_dir>/<basename>.gz.
    """
    dst_dir.mkdir(parents=True, exist_ok=True)

    if src.suffix == ".gz":
        dst = dst_dir / src.name
        if src.resolve() == dst.resolve():
            return dst
        shutil.copy2(src, dst)
        return dst

    # Not gz: gzip to dst_dir with .gz extension
    dst = dst_dir / (src.name + ".gz")
    with src.open("rb") as f_in, gzip.open(dst, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    return dst

def validate_and_process_fasta_files(data: Dict, fasta_files: List[str], outdir: Path) -> Dict:
    """
    Validate that all required FASTA files exist; place gzipped copies in outdir/fasta.
    Supports workbook entries that are either absolute/relative paths or basenames.
    """
    # Build indices for matching
    provided_paths: List[Path] = [Path(p).expanduser().resolve() for p in fasta_files]
    for p in provided_paths:
        if not p.exists():
            sys.exit(f"Error: FASTA file not found: {p}")

    by_full: Dict[str, Path] = {str(p): p for p in provided_paths}
    by_base: Dict[str, List[Path]] = defaultdict(list)
    for p in provided_paths:
        by_base[p.name].append(p)

    # Detect duplicate basenames (ambiguous)
    dupes = {b: lst for b, lst in by_base.items() if len(lst) > 1}
    if dupes:
        msg = "\n".join([f"- {b}:\n  " + "\n  ".join(map(str, lst)) for b, lst in dupes.items()])
        sys.exit(
            "Error: Duplicate basenames among provided FASTA files (ambiguous references in the workbook).\n"
            "Please disambiguate by using unique filenames or supply full paths in the workbook.\n" + msg
        )

    fasta_outdir = outdir / "fasta"

    for sra_id, segments in data.items():
        for segment, info in segments.items():
            ref = str(info.get('file') or "").strip()
            if not ref:
                logger.warning(f"Empty filename for {sra_id}/{segment}")
                continue

            # Try exact path match first, then basename match
            candidate = None
            ref_path = Path(ref)
            # If workbook provided a path, try to resolve relative to CWD
            if ref_path.exists():
                candidate = ref_path.resolve()
            elif ref in by_full:
                candidate = by_full[ref]
            else:
                base = Path(ref).name
                if base in by_base:
                    candidate = by_base[base][0]

            if candidate is None or not candidate.exists():
                sys.exit(f'Error: Required file "{ref}" (for {sra_id}/{segment}) was not found among the provided --fasta files.')

            staged = stage_gzip_to_outdir(candidate, fasta_outdir)
            data[sra_id][segment]['file'] = str(staged)

    return data

def export_data(data: Dict, submitter: str, outdir: Path) -> None:
    json_filename = outdir / f'{sanitize_string(submitter)}.json'
    try:
        with json_filename.open('w', encoding='utf-8') as f:
            json.dump(data, f, indent=4, ensure_ascii=False)
        logger.info(f'Data exported to {json_filename}')
    except Exception as e:
        logger.error(f"Failed to export JSON: {e}")
        sys.exit(f"Error: Could not write {json_filename}")

    csv_rows = []
    for sra_id, segments in data.items():
        for segment, info in segments.items():
            csv_rows.append({
                'sample': f'{sanitize_string(sra_id)}_{sanitize_string(segment)}',
                'assembly': info['file'],     # staged path in outdir/fasta
                'submitter': submitter
            })

    csv_filename = outdir / f'{sanitize_string(submitter)}.csv'
    try:
        with csv_filename.open('w', newline='', encoding='utf-8') as f:
            if csv_rows:
                writer = csv.DictWriter(f, fieldnames=csv_rows[0].keys())
                writer.writeheader()
                writer.writerows(csv_rows)
        logger.info(f'Data exported to {csv_filename}')
    except Exception as e:
        logger.error(f"Failed to export CSV: {e}")
        sys.exit(f"Error: Could not write {csv_filename}")

def main():
    parser = argparse.ArgumentParser(
        description="Process genomic assembly metadata and FASTA files for submission preparation",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--meta", type=str, required=True, help="Path to Excel metadata file")
    parser.add_argument("--fasta", type=str, required=True, nargs='+', help="FASTA file paths (space-separated)")
    parser.add_argument("--outdir", type=str, default=".", help="Directory to write outputs and staged FASTAs (default: current directory)")
    parser.add_argument('--version', action='version', version=f'prepare-entry.py v{VERSION}')
    args = parser.parse_args()

    outdir = Path(args.outdir).expanduser().resolve()
    ensure_outdir(outdir)

    startup_message = f"""
    prepare-entry.py v{VERSION}
    Output directory: {outdir}
    Processing genomic assembly data...
    """
    print(textwrap.dedent(startup_message), flush=True)

    if not os.path.exists(args.meta):
        sys.exit(f"Error: Metadata file not found: {args.meta}")

    for fasta_file in args.fasta:
        if not os.path.exists(fasta_file):
            sys.exit(f"Error: FASTA file not found: {fasta_file}")

    try:
        logger.info("Loading workbook data...")
        submitter, data = load_workbook_data(args.meta)

        logger.info(f"Found data for submitter: {submitter}")
        logger.info(f"Processing {len(data)} SRA entries...")

        logger.info("Validating and staging FASTA files to outdir...")
        data = validate_and_process_fasta_files(data, args.fasta, outdir)

        logger.info("Exporting processed data...")
        export_data(data, submitter, outdir)

        logger.info("Processing completed successfully!")

    except KeyboardInterrupt:
        logger.info("Processing interrupted by user")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()