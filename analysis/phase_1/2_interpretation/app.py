import streamlit as st
from pathlib import Path
from typing import List, Tuple
import numpy as np
import pandas as pd
import plotly.express as px

import io_ops
import stats

st.set_page_config(layout="wide")

# ============================================================================
# Load Data
# ============================================================================

df = None
df_filtered = None

st.header("Load Metric Files")
metric_paths = st.file_uploader(
    "Load metrics",
    type=["csv"],
    accept_multiple_files=True
)
if metric_paths:
    df = io_ops.load_stats_csv(metric_paths)

st.divider()

# ============================================================================
# Helper Functions
# ============================================================================

def _editor_with_checkboxes(
    df_in: pd.DataFrame,
    id_col: str,
    key: str,
    title: str
) -> list[str]:
    """Display data editor with checkboxes and return selected IDs."""
    if df_in is None or df_in.empty:
        st.warning(f"No data available for {title}.")
        return []
    
    if id_col not in df_in.columns:
        st.warning(f"Missing expected id column '{id_col}' in {title}.")
        return []

    view = df_in.copy()
    if "keep" not in view.columns:
        view["keep"] = True
    view = view[['keep'] + [c for c in df_in.columns if c != 'keep']]

    st.subheader(title)
    edited = st.data_editor(
        view,
        key=key,
        hide_index=True,
        use_container_width=True,
        column_config={
            "keep": st.column_config.CheckboxColumn("keep", default=True),
            id_col: st.column_config.TextColumn(id_col, disabled=True),
        }
    )
    
    keep_ids = edited.loc[edited["keep"] == True, id_col].astype(str).tolist()
    return keep_ids


def _filter_pairs_df(
    df_pairs: pd.DataFrame,
    keep_submitters: list[str] | None,
    keep_samples: list[str] | None
) -> pd.DataFrame:
    """Filter dataframe by submitters and samples."""
    out = df_pairs.copy()

    if keep_submitters and all(c in out.columns for c in ("seqA", "seqB")):
        out = out[
            out["seqA"].astype(str).isin(keep_submitters) &
            out["seqB"].astype(str).isin(keep_submitters)
        ]

    sample_cols = [
        c for c in ("sample", "id_alt", "sample_id", "id")
        if c in out.columns
    ]
    if keep_samples and sample_cols:
        sc = sample_cols[0]
        out = out[out[sc].astype(str).isin(keep_samples)]
    
    return out


def _sort_if_present(
    df_in: pd.DataFrame,
    col: str,
    ascending=True
) -> pd.DataFrame:
    """Sort dataframe by column if it exists."""
    if col in df_in.columns:
        return df_in.sort_values(col, ascending=ascending)
    return df_in

# ============================================================================
# Apply Filters
# ============================================================================

if df is not None and not df.empty:
    st.header("Apply Filters")

    # Compute row-level metrics
    df_rows = stats.row_stats(df)

    # Stage 1: Filter by valid_frac threshold
    st.markdown("### Filter 1: Minimum sample overlap")
    min_valid = st.number_input(
        "Minimum overlap (valid_frac) in sample assemblies between two submitters to be considered",
        min_value=0.0,
        max_value=1.0,
        step=0.01,
        value=0.80,
        format="%.2f",
        help="Rows with valid_frac > threshold are kept."
    )

    if "valid_frac" not in df_rows.columns:
        st.error("After row_stats, 'valid_frac' is missing. Please check stats.row_stats.")
        st.stop()

    df_stage1 = df_rows[df_rows["valid_frac"] > float(min_valid)].copy()
    st.caption(
        f"Kept {len(df_stage1)} / {len(df_rows)} rows "
        f"(valid_frac > {min_valid:.2f})."
    )

    # Build stats from Stage 1 filtered data
    df_submitter_stats = stats.submitter_stats(df_stage1)
    df_sample_stats = stats.sample_stats(df_stage1)

    # Stage 2: Checkbox selection for submitters and samples
    st.markdown("### Filter 2: Select submitters and samples")
    col1, col2 = st.columns([1, 1])

    with col1:
        if "submitter" not in df_submitter_stats.columns:
            st.error("submitter_stats must include a 'submitter' column.")
            st.stop()

        submitter_table = _sort_if_present(
            df_submitter_stats.drop_duplicates(),
            "match_frac_mean"
        )
        keep_submitters = _editor_with_checkboxes(
            df_in=submitter_table,
            id_col="submitter",
            key="editor_submitters",
            title="Select Submitters",
        )

    with col2:
        sample_id_col = (
            "sample" if "sample" in df_sample_stats.columns
            else df_sample_stats.columns[0]
        )
        sample_table = _sort_if_present(
            df_sample_stats.drop_duplicates(),
            "match_frac_mean"
        )
        keep_samples = _editor_with_checkboxes(
            df_in=sample_table,
            id_col=sample_id_col,
            key="editor_samples",
            title="Select Samples",
        )

    if not keep_submitters:
        st.info("No submitters selected — final dataset will be empty.")
    if not keep_samples:
        st.info("No samples selected — final dataset may be empty.")

    # Apply Stage 2 filtering
    df_filtered = _filter_pairs_df(
        df_stage1,
        keep_submitters=keep_submitters,
        keep_samples=keep_samples
    )

st.divider()

# ============================================================================
# Results: Plots & Statistics
# ============================================================================

if df_filtered is not None and not df_filtered.empty:
    st.header("Results")

    st.subheader("Global Stats (After filtering)")
    st.dataframe(stats.global_stats(df_filtered))

    st.subheader("Full Stats (After Filtering)")
    with st.expander("Expand to View", expanded=False):
        st.dataframe(df_filtered)

    coords_df, var_series, mat = stats.pcoa_classical(df_filtered)

    st.subheader("Plots")
    col1, col2 = st.columns([1, 1])
    
    with col1:
        fig_heatmap = stats.build_heatmap(mat)
        st.plotly_chart(fig_heatmap, use_container_width=True)

    with col2:
        fig_pcoa = stats.plot_pcoa(coords_df, var_series, "pcoa.html")
        st.plotly_chart(
            fig_pcoa,
            use_container_width=True,
            config={"displayModeBar": False}
        )

elif df is not None and df_filtered is not None and df_filtered.empty:
    st.warning(
        "No rows after Stage 2 filtering. "
        "Adjust your selections or lower the valid_frac threshold."
    )