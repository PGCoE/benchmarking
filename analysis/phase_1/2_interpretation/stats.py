from typing import List, Tuple, Iterable, Optional
from pathlib import Path
from collections import defaultdict

import numpy as np
import pandas as pd
import scipy.cluster.hierarchy as sch
from scipy.spatial.distance import squareform
import plotly.express as px
import plotly.graph_objects as go
import plotly.figure_factory as ff
from plotly.subplots import make_subplots

# ============================================================================
# Constants
# ============================================================================

NUMERICAL_COLS = [
    "match", "match_iupac", "ins", "del", "sub",
    "invalid", "valid", "valid_s1", "invalid_s1",
    "valid_s2", "invalid_s2", "del_terminal",
    "del_internal", "ins_terminal", "ins_internal"
]

STAT_COLS = [
    "match", "ins", "del", "sub", "valid"
]

# ============================================================================
# Row-Level Statistics
# ============================================================================

def row_stats(df):
    """
    Compute row-level statistics including valid fraction and metric fractions.
    Canonicalizes sequence pairs (seqA, seqB) for consistent ordering.
    """
    out = df.copy()

    # Convert numerical columns to numeric type
    for col in NUMERICAL_COLS:
        out[col] = pd.to_numeric(out[col], errors="coerce")

    # Canonicalize pairs: ensure seqA <= seqB
    a = out["seq1"].astype(str)
    b = out["seq2"].astype(str)
    seqA = a.where(a <= b, b)
    seqB = b.where(a <= b, a)
    out["seqA"] = seqA
    out["seqB"] = seqB

    # Fraction of valid sites based on total sites (valid + invalid)
    out["valid_frac"] = out["valid"] / (out["valid"] + out["invalid"])

    # Fraction of site classifications based on valid sites
    for col in STAT_COLS:
        if col in ["valid", "invalid"]:
            continue
        out[col + "_frac"] = out[col] / out["valid"]

    return out

# ============================================================================
# Submitter-Level Statistics
# ============================================================================

def submitter_stats(df: pd.DataFrame) -> pd.DataFrame:
    """
    Aggregate statistics per submitter across all their comparisons.
    Returns mean, standard deviation, and median for each metric.
    """
    stats_cols = [c + "_frac" for c in STAT_COLS]
    agg = defaultdict(lambda: defaultdict(list))

    # Collect values per submitter
    for rec in df.to_dict(orient="records"):
        for submitter in (rec["seqA"], rec["seqB"]):
            for metric in stats_cols:
                val = rec.get(metric)
                if pd.notna(val):
                    agg[submitter][metric].append(val)

    # Build output rows
    out_rows = []
    for submitter, metrics in agg.items():
        row = {"submitter": submitter}

        # Count total contributing samples
        sample_counts = [len(v) for v in metrics.values()]
        row["n_samples"] = int(np.nanmax(sample_counts)) if sample_counts else 0

        # Compute statistics for each metric
        for key, vals in metrics.items():
            arr = np.asarray(vals, dtype="float64")
            n = arr.size
            row[f"{key}_mean"] = float(np.nanmean(arr)) if n else np.nan
            row[f"{key}_sd"] = float(np.nanstd(arr, ddof=1)) if n > 1 else np.nan
            row[f"{key}_median"] = float(np.nanmedian(arr)) if n else np.nan
        
        out_rows.append(row)

    # Handle empty case
    if not out_rows:
        return pd.DataFrame(columns=["submitter", "n_samples"] + [
            f"{m}_frac_{suf}" for m in stats_cols for suf in ("mean", "sd", "median")
        ])

    df_out = pd.DataFrame(out_rows)

    # Establish stable column order
    ordered_cols = ["submitter", "n_samples"]
    for m in stats_cols:
        base = f"{m}_frac"
        ordered_cols += [f"{base}_mean", f"{base}_sd", f"{base}_median"]

    df_out = df_out[
        [c for c in ordered_cols if c in df_out.columns] +
        [c for c in df_out.columns if c not in ordered_cols]
    ]

    return df_out.sort_values("submitter").reset_index(drop=True)

# ============================================================================
# Sample-Level Statistics
# ============================================================================

def sample_stats(df: pd.DataFrame) -> pd.DataFrame:
    """
    Aggregate statistics per sample across all submitter comparisons.
    Returns mean, standard deviation, and median for each metric.
    """
    stats_cols = [c + "_frac" for c in STAT_COLS]
    agg = defaultdict(lambda: defaultdict(list))

    # Collect values per sample
    for rec in df.to_dict(orient="records"):
        sample = rec.get("sample")
        if not sample:
            continue
        for metric in stats_cols:
            val = rec.get(metric)
            if pd.notna(val):
                agg[sample][metric].append(val)

    # Build output rows
    out_rows = []
    for sample, metrics in agg.items():
        row = {"sample": sample}

        # Count total contributing submitters
        submitter_counts = [len(v) for v in metrics.values()]
        row["n_submitters"] = (
            int(np.nanmax(submitter_counts)) if submitter_counts else 0
        )

        # Compute statistics for each metric
        for key, vals in metrics.items():
            arr = np.asarray(vals, dtype="float64")
            n = arr.size
            row[f"{key}_mean"] = float(np.nanmean(arr)) if n else np.nan
            row[f"{key}_sd"] = float(np.nanstd(arr, ddof=1)) if n > 1 else np.nan
            row[f"{key}_median"] = float(np.nanmedian(arr)) if n else np.nan
        
        out_rows.append(row)

    # Handle empty case
    if not out_rows:
        return pd.DataFrame(columns=["sample", "n_submitters"] + [
            f"{m}_frac_{suf}" for m in stats_cols for suf in ("mean", "sd", "median")
        ])

    df_out = pd.DataFrame(out_rows)

    # Establish stable column order
    ordered_cols = ["sample", "n_submitters"]
    for m in stats_cols:
        base = f"{m}_frac"
        ordered_cols += [f"{base}_mean", f"{base}_sd", f"{base}_median"]

    df_out = df_out[
        [c for c in ordered_cols if c in df_out.columns] +
        [c for c in df_out.columns if c not in ordered_cols]
    ]

    return df_out.sort_values("sample").reset_index(drop=True)

# ============================================================================
# Global Statistics
# ============================================================================

def global_stats(df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute global statistics across all pairwise comparisons.
    Returns min, max, mean, standard deviation, and median for each metric.
    """
    stats_cols = [c + "_frac" for c in STAT_COLS]
    agg = defaultdict(list)

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

    return df_out

# ============================================================================
# Distance Matrix Operations
# ============================================================================

def to_distance_matrix(agg: pd.DataFrame) -> pd.DataFrame:
    """
    Build a symmetric square matrix from aggregated pair dissimilarities.
    Diagonal set to 0.0. Missing pairs will have NaN values.
    """
    samples = sorted(set(agg["seqA"]).union(set(agg["seqB"])))
    mat = pd.DataFrame(index=samples, columns=samples, dtype=float)
    
    # Set diagonal to zero
    np.fill_diagonal(mat.values, 0.0)
    
    # Fill both [A,B] and [B,A] for symmetry
    for _, row in agg.iterrows():
        a, b, d = row["seqA"], row["seqB"], row["mean_dissimilarity_frac"]
        mat.at[a, b] = d
        mat.at[b, a] = d
    
    return mat


def mat_drop_na(mat: pd.DataFrame) -> pd.DataFrame:
    """
    Iteratively drop rows/columns with the most NaNs until a complete
    submatrix remains. Maintains square structure throughout.
    """
    current = mat.copy()

    while True:
        # Count NaNs per row and column
        row_na = (
            current.isna().sum(axis=1) if current.shape[0]
            else pd.Series(dtype=int)
        )
        col_na = (
            current.isna().sum(axis=0) if current.shape[1]
            else pd.Series(dtype=int)
        )
        r_max = int(row_na.max()) if not row_na.empty else 0
        c_max = int(col_na.max()) if not col_na.empty else 0

        # Exit if no NaNs remain
        if r_max == 0 and c_max == 0:
            break

        # Drop the worst offender (row ties favor rows)
        if r_max >= c_max:
            drop_label = row_na.idxmax()
            current = current.drop(index=drop_label, errors="ignore")
        else:
            drop_label = col_na.idxmax()
            current = current.drop(columns=drop_label, errors="ignore")

        # Re-square to intersection of remaining labels
        keep = current.index.intersection(current.columns)
        current = current.loc[keep, keep]

        if current.shape[0] == 0:
            raise ValueError(
                "No fully observed submatrix exists (too many missing pairs)."
            )

    return current

# ============================================================================
# Principal Coordinates Analysis (PCoA)
# ============================================================================

def pcoa_classical(
    df: pd.DataFrame,
    k: int = 2
) -> Tuple[pd.DataFrame, pd.Series, pd.DataFrame]:
    """
    Classical multidimensional scaling (PCoA) from a metric distance matrix.
    
    Returns:
        coords_df: DataFrame with PCo coordinates (n x k)
        var_series: Series with percent variance explained per axis
        distance_matrix: The processed distance matrix
    """
    df_subset = df[["seqA", "seqB", "match_frac"]]
    df_subset["dissimilarity_frac"] = 1 - df["match_frac"]

    # Calculate average dissimilarity per pair
    agg = (
        df_subset.groupby(["seqA", "seqB"], as_index=False)
        .agg(mean_dissimilarity_frac=("dissimilarity_frac", "mean"))
    )

    distance_matrix = to_distance_matrix(agg)
    distance_matrix = mat_drop_na(distance_matrix)

    D = distance_matrix.values.astype(float)
    
    # Validate matrix properties
    if not np.allclose(D, D.T, equal_nan=False):
        raise ValueError("Distance matrix is not symmetric.")
    if np.any(np.diag(D) != 0):
        raise ValueError("Distance matrix diagonal must be zero.")

    # Double-center the squared distances
    n = D.shape[0]
    J = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * J @ (D ** 2) @ J

    # Eigen decomposition
    eigvals, eigvecs = np.linalg.eigh(B)
    
    # Sort eigenvalues/vectors in descending order
    idx = np.argsort(eigvals)[::-1]
    eigvals = eigvals[idx]
    eigvecs = eigvecs[:, idx]

    # Keep only positive eigenvalues
    pos = eigvals > 1e-12
    eigvals_pos = eigvals[pos]
    eigvecs_pos = eigvecs[:, pos]

    if eigvals_pos.size == 0:
        raise ValueError("No positive eigenvalues; PCoA failed.")

    # Calculate coordinates
    L = np.sqrt(eigvals_pos[:k])
    coords = eigvecs_pos[:, :k] * L

    # Variance explained
    var_explained = eigvals_pos / eigvals_pos.sum()
    
    coords_df = pd.DataFrame(
        coords,
        index=distance_matrix.index,
        columns=[f"PCo{c+1}" for c in range(coords.shape[1])]
    )
    var_series = pd.Series(
        var_explained[:coords.shape[1]],
        index=coords_df.columns
    )
    
    return coords_df, var_series, distance_matrix

# ============================================================================
# Visualization Functions
# ============================================================================

def plot_pcoa(
    coords: pd.DataFrame,
    var: pd.Series,
    out_html: Path
) -> go.Figure:
    """
    Create an interactive Plotly PCoA scatterplot.
    
    Args:
        coords: DataFrame with PCo coordinates (index = sample IDs)
        var: Series with percent variance explained per axis
        out_html: Output HTML file path
    
    Returns:
        Plotly Figure object
    """
    # Prepare plotting DataFrame
    df_plot = coords.copy()
    df_plot["sample"] = df_plot.index

    # Axis labels with variance explained
    xcol = coords.columns[0]
    ycol = coords.columns[1] if coords.shape[1] > 1 else None

    xlab = f"{xcol} ({var.iloc[0]*100:.1f}% var)"
    if ycol:
        ylab = f"{ycol} ({var.iloc[1]*100:.1f}% var)"
    else:
        ylab = ycol or "PCo2"

    # Build Plotly figure
    fig = px.scatter(
        df_plot,
        x=xcol,
        y=ycol,
        text="sample",
        hover_name="sample",
        hover_data=df_plot.columns.tolist(),
        title="PCoA (classical MDS) on pairwise dissimilarities",
    )

    # Customize appearance
    fig.update_traces(textposition="top center", marker=dict(size=10))
    fig.update_layout(
        xaxis_title=xlab,
        yaxis_title=ylab,
        height=650,
        width=900,
        margin=dict(l=40, r=40, t=60, b=60)
    )

    # Save to HTML
    fig.write_html(str(out_html), include_plotlyjs="cdn")

    return fig


def build_heatmap(
    mat: pd.DataFrame,
    method: str = "average",
    title: str = "Hierarchical Heatmap"
) -> go.Figure:
    """
    Build a hierarchical clustered heatmap with top dendrogram.
    
    Args:
        mat: Square distance matrix with matching index and columns
        method: Linkage method for hierarchical clustering
        title: Plot title
    
    Returns:
        Plotly Figure object
    """
    # Validate input
    assert mat.shape[0] == mat.shape[1], f"Square matrix required, got {mat.shape}"
    assert mat.index.equals(mat.columns), (
        "Index and columns must be identical and aligned."
    )
    labels = mat.index.astype(str).tolist()

    # Convert to numeric array
    A = mat.apply(pd.to_numeric, errors="coerce").to_numpy(
        dtype=np.float64, copy=False
    )
    
    # Ensure zero diagonal
    np.fill_diagonal(A, 0.0)
    
    # Replace non-finite values
    if not np.isfinite(A).all():
        finite = A[np.isfinite(A)]
        repl = (np.nanmax(finite) if finite.size else 1.0) * 10.0
        A = np.where(np.isfinite(A), A, repl).astype(np.float64, copy=False)

    # Perform hierarchical clustering
    condensed = squareform(A, checks=False).astype(np.float64, copy=False)
    Z = sch.linkage(condensed, method=method)
    Z = np.asarray(Z, dtype=np.float64, order="C")
    order = sch.leaves_list(Z)

    # Reorder matrix and labels
    A_ord = A[np.ix_(order, order)]
    labels_ord = [labels[i] for i in order]

    # Build top dendrogram
    dendro_top = ff.create_dendrogram(
        A_ord,
        orientation="top",
        labels=labels_ord,
        distfun=lambda _x: squareform(A_ord, checks=False),
        linkagefun=lambda d: sch.linkage(d, method=method),
    )
    for t in dendro_top.data:
        t.showlegend = False

    # Create subplots: dendrogram + heatmap
    fig = make_subplots(
        rows=2, cols=1,
        row_heights=[0.22, 0.78],
        specs=[[{"type": "xy"}], [{"type": "heatmap"}]],
        vertical_spacing=0.01,
    )

    # Add dendrogram traces
    for tr in dendro_top.data:
        fig.add_trace(tr, row=1, col=1)

    # Add heatmap
    heat = go.Heatmap(
        z=A_ord,
        x=labels_ord,
        y=labels_ord,
        colorscale="Viridis",
        colorbar=dict(title="Dissimilarity")
    )
    fig.add_trace(heat, row=2, col=1)

    # Configure axes
    # Dendrogram facing downward
    fig.update_yaxes(autorange="reversed", row=1, col=1)
    fig.update_xaxes(autorange="reversed", row=1, col=1)
    fig.update_xaxes(showticklabels=False, row=1, col=1)
    fig.update_yaxes(visible=False, row=1, col=1)

    # Heatmap: y-axis reversed so first row is at top
    fig.update_xaxes(side="bottom", row=2, col=1)
    fig.update_yaxes(autorange="reversed", row=2, col=1)

    fig.update_layout(
        width=1000,
        height=850,
        title=title,
        margin=dict(l=10, r=10, t=40, b=10),
    )
    
    return fig