import pandas as pd
import streamlit as st
from collections import defaultdict
import numpy as np

METRICS_COLS = [
    "match", "ins", "del", "sub",
    "invalid", "valid", "valid_s1", "invalid_s1",
    "valid_s2", "invalid_s2", "del_terminal",
    "del_internal", "ins_terminal", "ins_internal",
    "invalid_char_s1", "invalid_char_s2"
]

STAT_COLS = [
    "match", "ins", "del", "sub", "valid"
]

def add_attr(df_metrics, df_attr):
    # --- merge for seq1 ---
    out = df_metrics.merge(
        df_attr,
        left_on=["sample", "source_s1"],
        right_on=["sample", "source"],
        how="left",
        suffixes=("", "_1")
    )

    # rename merged columns from first merge
    meta_cols = df_attr.columns.difference(["sample", "source"])
    out = out.rename(columns={c: f"{c}_s1" for c in meta_cols})

    # --- merge for seq2 ---
    out = out.merge(
        df_attr,
        left_on=["sample", "source_s2"],
        right_on=["sample", "source"],
        how="left"
    )

    # rename second merge columns
    out = out.rename(columns={c: f"{c}_s2" for c in meta_cols})

    # optional: drop duplicate source columns
    out = out.drop(columns=["source_x", "source_y"], errors="ignore")

    return out

def row_stats(df, ignore_termini: bool = False):
    out = df.copy()
    # Convert numerical columns to numeric type
    for col in METRICS_COLS:
        out[col] = pd.to_numeric(out[col], errors="coerce")

    if ignore_termini:
        term = out['ins_terminal'] + out['del_terminal']

        out['ins'] = out['ins'] - out['ins_terminal']
        out['del'] = out['del'] - out['del_terminal']

        out['valid']   = out['valid'] - term
        out['invalid'] = out['invalid'] + term

    # Fraction of valid sites based on total sites (valid + invalid)
    out["valid_frac"] = out["valid"] / (out["valid"] + out["invalid"])

    # Fraction of site classifications based on valid sites
    for col in STAT_COLS:
        if col in ["valid", "invalid"]:
            continue
        out[col + "_frac"] = out[col] / out["valid"]

    return out

def apply_grouping(col_in, col_out):
    df = st.session_state.df
    col_1, col_2 = f"{col_in}_s1", f"{col_in}_s2"
    if col_1 not in df.columns or col_2 not in df.columns:
        col_1 = col_2 = col_in

    a = df[col_1].astype(str)
    b = df[col_2].astype(str)
    df[col_out + '_a'] = a.where(a <= b, b)
    df[col_out + '_b'] = b.where(a <= b, a)

    st.session_state.df = df


# ============================================================================
# Submitter-Level Statistics
# ============================================================================

def group_stats() -> pd.DataFrame:
    df = st.session_state.df_final

    stats_cols = [c + "_frac" for c in METRICS_COLS]

    agg = defaultdict(lambda: defaultdict(list))
    samples_by_group = defaultdict(set)

    for rec in df.to_dict(orient="records"):
        sample = rec.get("sample")

        for member in (rec["gcol_a"], rec["gcol_b"]):
            if pd.notna(sample):
                samples_by_group[member].add(sample)

            for metric in stats_cols:
                val = rec.get(metric)
                if pd.notna(val):
                    agg[member][metric].append(val)

    out_rows = []
    for member, metrics in agg.items():
        row = {
            "member": member,
            "n_samples": len(samples_by_group.get(member, set())),
        }

        for key, vals in metrics.items():
            arr = np.asarray(vals, dtype="float64")
            n = arr.size
            row[f"{key}_mean"] = float(np.nanmean(arr)) if n else np.nan
            row[f"{key}_sd"] = float(np.nanstd(arr, ddof=1)) if n > 1 else np.nan
            row[f"{key}_median"] = float(np.nanmedian(arr)) if n else np.nan

        out_rows.append(row)

    if not out_rows:
        return pd.DataFrame(columns=["member", "n_samples"] + [
            f"{m}_{suf}" for m in stats_cols for suf in ("mean", "sd", "median")
        ])

    df_out = pd.DataFrame(out_rows)

    ordered_cols = ["member", "n_samples"]
    for m in stats_cols:
        ordered_cols += [f"{m}_mean", f"{m}_sd", f"{m}_median"]

    df_out = df_out[
        [c for c in ordered_cols if c in df_out.columns] +
        [c for c in df_out.columns if c not in ordered_cols]
    ]

    print(df_out.sort_values("member").reset_index(drop=True))


def to_distance_matrix() -> pd.DataFrame:
    df = st.session_state.df_final[["gcol_a", "gcol_b", "match_frac"]]
    df["dissimilarity_frac"] = 1 - df["match_frac"]

    # Calculate average dissimilarity per pair
    agg = (
        df.groupby(["gcol_a", "gcol_b"], as_index=False)
        .agg(mean_dissimilarity_frac=("dissimilarity_frac", "mean"))
    )

    samples = sorted(set(agg["gcol_a"]).union(set(agg["gcol_b"])))
    mat = pd.DataFrame(index=samples, columns=samples, dtype=float)
    
    # Set diagonal to zero
    np.fill_diagonal(mat.values, 0.0)
    
    # Fill both [A,B] and [B,A] for symmetry
    for _, row in agg.iterrows():
        a, b, d = row["gcol_a"], row["gcol_b"], row["mean_dissimilarity_frac"]
        mat.at[a, b] = d
        mat.at[b, a] = d

    D = mat.values.astype(float)

    if not np.allclose(D, D.T, equal_nan=False):
        raise ValueError("Distance matrix is not symmetric.")

    st.session_state.matrix_values = D
    st.session_state.matrix_labels = samples


def global_stats() -> pd.DataFrame:
    """
    Compute global statistics across all pairwise comparisons.
    Returns min, max, mean, standard deviation, and median for each metric.
    """
    stats_cols = [c + "_frac" for c in STAT_COLS]
    agg = defaultdict(list)

    df = st.session_state.df_final

    # Collect all values per metric
    for rec in df.to_dict(orient="records"):
        for metric in stats_cols:
            val = rec.get(metric)
            if pd.notna(val):
                agg[metric].append(val)

    # Compute summary statistics per metric
    summary = {}
    for metric, vals in agg.items():
        arr = np.asarray(vals, dtype="float64")
        if arr.size == 0:
            summary[metric] = dict(
                min=np.nan, max=np.nan, mean=np.nan, sd=np.nan, median=np.nan
            )
        else:
            summary[metric] = dict(
                min=float(np.nanmin(arr)),
                max=float(np.nanmax(arr)),
                mean=float(np.nanmean(arr)),
                sd=float(np.nanstd(arr, ddof=1)) if arr.size > 1 else np.nan,
                median=float(np.nanmedian(arr)),
            )

    # Convert to DataFrame with statistic rows
    df_out = pd.DataFrame(summary).rename_axis(index="statistic")
    df_out = df_out.loc[["min", "max", "mean", "sd", "median"]]

    st.session_state.df_global = df_out

