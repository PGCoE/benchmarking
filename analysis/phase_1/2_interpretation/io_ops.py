from pathlib import Path
from typing import List
import pandas as pd

# ============================================================================
# Constants
# ============================================================================

INPUT_STATS_COLS = [
    "match", "match_iupac", "ins", "del", "sub",
    "invalid", "valid", "valid_s1", "invalid_s1",
    "valid_s2", "invalid_s2", "del_terminal",
    "del_internal", "ins_terminal", "ins_internal",
    "seq1", "seq2", "invalid_char_s1", 
    "invalid_char_s2", "sample"
]

# ============================================================================
# Data Loading Functions
# ============================================================================

def load_stats_csv(csv_paths: List[Path]) -> pd.DataFrame:
    """
    Load and concatenate multiple CSV files containing alignment statistics.
    
    Args:
        csv_paths: List of paths to CSV files
    
    Returns:
        Concatenated DataFrame with all statistics
    
    Raises:
        ValueError: If required columns are missing or no data is found
    """
    dfs = []
    
    for p in csv_paths:
        df = pd.read_csv(p)

        # Validate required columns are present
        missing = set(INPUT_STATS_COLS) - set(df.columns)
        if missing:
            raise ValueError(
                f"{p}: missing required columns: {sorted(missing)}"
            )

        dfs.append(df)

    if not dfs:
        raise ValueError("No input rows found.")

    return pd.concat(dfs, ignore_index=True)