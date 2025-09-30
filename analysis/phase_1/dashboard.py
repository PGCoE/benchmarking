#!/usr/bin/env python3

# dashboard.py
# Author: Jared Johnson, jared.johnson@doh.wa.gov

"""
Streamlit app for interpreting viral benchmarking metrics
"""

import io
from collections import defaultdict

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from sklearn.decomposition import PCA
import streamlit as st


# ----------------------------
# Helpers
# ----------------------------
NUM_KEYS_SUBMITTER = (
    "valid", "match", "match_iupac", "ins", "del", "del_internal", "del_terminal",
    "ins_internal", "ins_terminal", "sub"
)

def num(row, key, default=0):
    """Safe numeric conversion with default 0 for NA/empty."""
    v = pd.to_numeric(row.get(key, default), errors="coerce")
    return 0 if (v is None or pd.isna(v)) else v

def pct_of(part, whole):
    return (100.0 * part / whole) if whole > 0 else 0.0

def effective_match_sub(match, match_iupac, sub, use_iupac=False):
    """Return (effective_match, effective_sub) given IUPAC policy."""
    if use_iupac:
        return match + match_iupac, sub
    return match, sub + match_iupac

def calculate_error_rate(match, valid, match_iupac=0, ins=0, delet=0, sub=0, 
                         ins_terminal=0, del_terminal=0, use_iupac=False, 
                         include_terminal=True):
    """Calculate error rate with optional terminal indel exclusion."""
    if not valid or pd.isna(valid) or valid == 0:
        return 0.0
    
    eff_match, eff_sub = effective_match_sub(match, match_iupac, sub, use_iupac)
    
    # Calculate errors
    if include_terminal:
        errors = valid - eff_match
    else:
        # Exclude terminal indels from error calculation
        errors = (eff_sub + (ins - ins_terminal) + (delet - del_terminal))
    
    return max(0.0, errors / valid)

def load_csv_files(uploaded_files):
    """Load and combine multiple CSV files into list of dict rows with stripped column names."""
    rows = []
    for f in uploaded_files:
        df = pd.read_csv(io.StringIO(f.read().decode("utf-8")))
        f.seek(0)
        df.columns = df.columns.str.strip()
        rows.extend(df.to_dict("records"))
    return rows


# ----------------------------
# Matrix building + preprocessing
# ----------------------------
def build_pca_matrix(rows, use_iupac=False, include_terminal=True):
    """Build comparison matrix (rows=comparisons; cols=samples) of error rates."""
    comparisons = defaultdict(dict)
    samples = set()

    for r in rows:
        sample = r.get("sample")
        seq1, seq2 = str(r.get("seq1", "")), str(r.get("seq2", ""))
        if not (sample and seq1 and seq2):
            continue

        match = num(r, "match")
        valid = num(r, "valid")
        match_iupac = num(r, "match_iupac")
        ins = num(r, "ins")
        delet = num(r, "del")
        sub = num(r, "sub")
        ins_terminal = num(r, "ins_terminal")
        del_terminal = num(r, "del_terminal")

        er = calculate_error_rate(match, valid, match_iupac, ins, delet, sub,
                                 ins_terminal, del_terminal, use_iupac, include_terminal)
        comparisons[(seq1, seq2)][sample] = er
        samples.add(sample)

    row_keys = sorted(comparisons.keys())
    col_keys = sorted(samples)
    if not row_keys or not col_keys:
        return np.empty((0, 0)), row_keys, col_keys

    M = np.full((len(row_keys), len(col_keys)), np.nan)
    for i, rk in enumerate(row_keys):
        for j, ck in enumerate(col_keys):
            if ck in comparisons[rk]:
                M[i, j] = comparisons[rk][ck]
    return M, row_keys, col_keys

def impute_and_clean_matrix(matrix, col_keys):
    """Impute NaNs with column means; drop constant columns."""
    if matrix.size == 0:
        return matrix, col_keys

    # Impute NaNs with column means (or 0 if column all-NaN)
    col_means = np.nanmean(matrix, axis=0)
    col_means = np.where(np.isfinite(col_means), col_means, 0.0)
    inds = np.where(~np.isfinite(matrix))
    matrix[inds] = np.take(col_means, inds[1])

    # Remove near-constant columns
    std = np.nanstd(matrix, axis=0)
    keep = std > 1e-10
    if keep.any():
        matrix = matrix[:, keep]
        col_keys = [c for k, c in enumerate(col_keys) if keep[k]]
    else:
        matrix = np.empty((matrix.shape[0], 0))
        col_keys = []
    return matrix, col_keys


# ----------------------------
# Stats tables
# ----------------------------
def calculate_submitter_stats(rows, use_iupac=False, include_terminal=True):
    """
    Aggregate per submitter (each unique sequence id acts as 'submitter').
    Keep original metric semantics and directional indel swapping for seq1 vs seq2.
    """
    totals = defaultdict(lambda: dict(rows=0, **{k: 0 for k in NUM_KEYS_SUBMITTER}))

    for r in rows:
        seq1, seq2 = r.get("seq1"), r.get("seq2")
        if not (seq1 and seq2):
            continue

        v = {k: num(r, k) for k in NUM_KEYS_SUBMITTER}

        sample = r.get('sample')

        # Update for seq2 (as written/original)
        t2 = totals[seq2]
        t2.setdefault('samples', set()).add(sample)
        t2["rows"] += 1
        for k in NUM_KEYS_SUBMITTER:
            t2[k] += v[k]

        # For seq1, swap ins<->del (and their internal/terminal breakdowns)
        t1 = totals[seq1]
        t1.setdefault('samples', set()).add(sample)
        t1["rows"] += 1
        t1["valid"]        += v["valid"]
        t1["match"]        += v["match"]
        t1["match_iupac"]  += v["match_iupac"]
        t1["ins"]          += v["del"]
        t1["ins_internal"] += v["del_internal"]
        t1["ins_terminal"] += v["del_terminal"]
        t1["del"]          += v["ins"]
        t1["del_internal"] += v["ins_internal"]
        t1["del_terminal"] += v["ins_terminal"]
        t1["sub"]          += v["sub"]

    records = []
    for submitter, t in totals.items():
        eff_match, eff_sub = effective_match_sub(t["match"], t["match_iupac"], t["sub"], use_iupac)

        if include_terminal:
            valid = t["valid"]
            ins = t["ins"];      _del = t["del"]
            ins_inter = t["ins_internal"]; del_inter = t["del_internal"]
            ins_term  = t["ins_terminal"]; del_term  = t["del_terminal"]
        else:
            # Exclude terminal indels from denominators and from component percentages
            term_sum = t["ins_terminal"] + t["del_terminal"]
            valid = max(0, t["valid"] - term_sum)
            ins = t["ins_internal"];      _del = t["del_internal"]
            ins_inter = t["ins_internal"]; del_inter = t["del_internal"]
            ins_term = 0;                 del_term = 0

        rec = {
            "submitter": submitter,
            "samples": len(t["samples"]),
            "rows": t["rows"],
            "valid": valid,
            "pct_match": pct_of(eff_match, valid),
            "pct_ins": pct_of(ins, valid),
            "pct_ins_internal": pct_of(ins_inter, valid),
            "pct_ins_terminal": pct_of(ins_term, valid),
            "pct_del": pct_of(_del, valid),
            "pct_del_internal": pct_of(del_inter, valid),
            "pct_del_terminal": pct_of(del_term, valid),
            "pct_sub": pct_of(eff_sub, valid),

            # raw counts preserved (unchanged totals)
            "ins_internal": t["ins_internal"],
            "ins_terminal": t["ins_terminal"],
            "del_internal": t["del_internal"],
            "del_terminal": t["del_terminal"],
            "include": True,
        }
        records.append(rec)

    df = pd.DataFrame.from_records(records)
    return df.sort_values(["pct_match", "rows"], ascending=[False, False]) if not df.empty else df


def calculate_sample_stats(rows, use_iupac=False, include_terminal=True):
    """Aggregate per-sample gap (ins+del), match, sub; preserve original math."""
    totals = defaultdict(lambda: {
        "rows": 0, "valid": 0, "match": 0, "match_iupac": 0,
        "gap": 0, "gap_internal": 0, "gap_terminal": 0, "sub": 0
    })

    for r in rows:
        sample = r.get("sample")
        if not sample:
            continue

        v_valid = num(r, "valid")
        v_match = num(r, "match")
        v_miup  = num(r, "match_iupac")
        v_ins   = num(r, "ins")
        v_del   = num(r, "del")
        v_iint  = num(r, "ins_internal")
        v_iterm = num(r, "ins_terminal")
        v_dint  = num(r, "del_internal")
        v_dterm = num(r, "del_terminal")
        v_sub   = num(r, "sub")

        if include_terminal:
            valid = v_valid
            ins = v_ins; _del = v_del
            gap_internal = v_iint + v_dint
            gap_terminal = v_iterm + v_dterm
        else:
            term_sum = v_iterm + v_dterm
            valid = max(0, v_valid - term_sum)
            ins = v_iint; _del = v_dint
            gap_internal = v_iint + v_dint
            gap_terminal = 0

        t = totals[sample]
        t["rows"] += 1
        t["valid"] += valid
        t["match"] += v_match
        t["match_iupac"] += v_miup
        t["gap"] += (ins + _del)
        t["gap_internal"] += gap_internal
        t["gap_terminal"] += gap_terminal
        t["sub"] += v_sub

    records = []
    for sample, t in totals.items():
        eff_match, eff_sub = effective_match_sub(t["match"], t["match_iupac"], t["sub"], use_iupac)
        records.append({
            "sample": sample,
            "rows": t["rows"],
            "valid": t["valid"],
            "pct_match": pct_of(eff_match, t["valid"]),
            "pct_gap": pct_of(t["gap"], t["valid"]),
            "pct_gap_internal": pct_of(t["gap_internal"], t["valid"]),
            "pct_gap_terminal": pct_of(t["gap_terminal"], t["valid"]),
            "pct_sub": pct_of(eff_sub, t["valid"]),
            "include": True,
        })

    df = pd.DataFrame.from_records(records)
    return df.sort_values(["pct_match", "rows"], ascending=[False, False]) if not df.empty else df



# ----------------------------
# Filtering
# ----------------------------
def filter_by_coverage(rows, selected_submitters, min_fraction):
    """Keep samples that have comparisons covering at least min_fraction of selected submitters."""
    if not selected_submitters:
        return []

    per_sample_submitters = defaultdict(set)
    sel = set(selected_submitters)
    for r in rows:
        s = r.get("sample")
        a, b = r.get("seq1"), r.get("seq2")
        if s and a in sel and b in sel:
            per_sample_submitters[s].update([a, b])

    n = len(sel)
    return [s for s, subs in per_sample_submitters.items() if len(subs) / n >= min_fraction]


# ----------------------------
# PCA + Plot
# ----------------------------
def create_pca_plot(matrix, row_keys, show_spokes=True, show_centroids=True, point_size=8):
    n_components = min(2, matrix.shape[1]) if matrix.size else 2
    pca = PCA(n_components=n_components)
    coords = pca.fit_transform(matrix)

    if coords.shape[1] == 1:  # pad PC2 if degenerate
        coords = np.hstack([coords, np.zeros((coords.shape[0], 1))])

    # group points by submitter (both seq1 and seq2 receive each comparison point)
    group_pts = defaultdict(list)
    for i, (s1, s2) in enumerate(row_keys):
        xy = coords[i, :2].tolist()
        group_pts[s1].append(xy)
        group_pts[s2].append(xy)

    # centroids
    centroids = {g: np.mean(pts, axis=0) for g, pts in group_pts.items() if pts}

    fig = go.Figure()
    palette = px.colors.qualitative.Set1
    color = {g: palette[i % len(palette)] for i, g in enumerate(sorted(group_pts))}

    # spokes (point -> centroid)
    if show_spokes:
        for g, pts in group_pts.items():
            if g not in centroids:
                continue
            cx, cy = centroids[g]
            xs, ys = [], []
            for x, y in pts:
                xs += [x, cx, None]
                ys += [y, cy, None]
            fig.add_trace(go.Scatter(x=xs, y=ys, mode="lines",
                                     line=dict(color=color[g], width=1),
                                     opacity=0.3, name=f"{g} spokes",
                                     showlegend=False, hoverinfo="skip"))

    # points
    for g, pts in group_pts.items():
        xs = [p[0] for p in pts]
        ys = [p[1] for p in pts]
        fig.add_trace(go.Scatter(x=xs, y=ys, mode="markers",
                                 marker=dict(color=color[g], size=point_size),
                                 opacity=0.7, name=g, showlegend=True))

    # centroids
    if show_centroids:
        for g, (cx, cy) in centroids.items():
            fig.add_trace(go.Scatter(x=[cx], y=[cy], mode="markers",
                                     marker=dict(color=color[g], size=18, symbol="circle",
                                                 line=dict(color="black", width=2)),
                                     name=f"{g} centroid", showlegend=True))

    varexp = pca.explained_variance_ratio_
    pc1 = round(100 * varexp[0], 2) if len(varexp) > 0 else 0.0
    pc2 = round(100 * varexp[1], 2) if len(varexp) > 1 else 0.0

    fig.update_layout(
        title="PCA of Submitter Comparisons",
        xaxis_title=f"PC1 ({pc1}%)",
        yaxis_title=f"PC2 ({pc2}%)",
        template="simple_white",
        width=900, height=650
    )
    return fig, pca

def calculate_variance_contributions(pca, col_keys):
    comps = pca.components_.T  # samples x PCs
    varexp = pca.explained_variance_ratio_
    if comps.size == 0:
        return pd.DataFrame(columns=["sample", "contribution", "percentage"])

    k = min(2, comps.shape[1])
    weights = varexp[:k]
    contrib = (comps[:, :k] ** 2) * weights
    total = contrib.sum(axis=1)
    total_sum = total.sum()

    df = pd.DataFrame({"sample": col_keys, "contribution": total})
    df["percentage"] = df["contribution"] / total_sum * 100.0 if total_sum > 0 else 0.0
    return df.sort_values("percentage", ascending=False)


# ----------------------------
# Streamlit UI
# ----------------------------
def main():
    st.set_page_config(page_title="Submitter PCA", layout="wide")
    st.title("Submitter Comparison PCA Analysis")
    st.caption("Upload CSVs to analyze submitter comparisons using Principal Component Analysis")

    with st.sidebar:
        st.header("Controls")
        files = st.file_uploader(
            "Upload CSV files",
            type=["csv"],
            accept_multiple_files=True,
            help="Files should contain: sample, seq1, seq2, match, valid columns"
        )
        use_iupac = st.checkbox(
            "Count IUPAC matches as true matches",
            value=False,
            help="Include match_iupac in match counts instead of substitutions"
        )
        include_terminal = st.checkbox(
            "Include terminal indels",
            value=True,
            help="Include terminal insertions and deletions in error calculations"
        )
        min_fraction = st.slider(
            "Minimum fraction of submitters per sample",
            min_value=0.0, max_value=1.0, value=1.0, step=0.01,
            help="Only keep samples covered by this fraction of submitters"
        )
        st.subheader("Display Options")
        show_spokes = st.checkbox("Show spokes to centroids", value=True)
        show_centroids = st.checkbox("Show centroids", value=True)
        point_size = st.slider("Point size", 4, 16, 8)

    if not files:
        st.info("Upload CSV files to begin analysis")
        return

    try:
        rows = load_csv_files(files)
        st.success(f"Loaded {len(rows)} rows from {len(files)} files")
    except Exception as e:
        st.error(f"Error loading files: {e}")
        return

    if not rows:
        st.error("No data found in uploaded files")
        return

    # Initial stats
    submitter_stats = calculate_submitter_stats(rows, use_iupac, include_terminal)
    sample_stats = calculate_sample_stats(rows, use_iupac, include_terminal)

    c1, c2 = st.columns(2)
    with c1:
        st.subheader("Submitter Statistics")
        show_cols = ["include", "submitter", "samples", "rows", "pct_match", "pct_ins", "pct_del", "pct_sub"]
        edited_submitters = st.data_editor(
            submitter_stats[show_cols] if not submitter_stats.empty else pd.DataFrame(columns=show_cols),
            use_container_width=True,
            hide_index=True
        )
        selected_submitters = set(edited_submitters.loc[edited_submitters.get("include", False) == True, "submitter"].tolist())

    with c2:
        st.subheader("Sample Statistics")
        show_cols = ["include", "sample", "rows", "pct_match", "pct_gap", "pct_sub"]
        edited_samples = st.data_editor(
            sample_stats[show_cols] if not sample_stats.empty else pd.DataFrame(columns=show_cols),
            use_container_width=True,
            hide_index=True
        )
        selected_samples = set(edited_samples.loc[edited_samples.get("include", False) == True, "sample"].tolist())

    # Filter rows by selections
    filtered_rows = [
        r for r in rows
        if (r.get("sample") in selected_samples
            and r.get("seq1") in selected_submitters
            and r.get("seq2") in selected_submitters)
    ]

    # Coverage filter
    kept_samples = filter_by_coverage(filtered_rows, selected_submitters, min_fraction)
    final_rows = [r for r in filtered_rows if r.get("sample") in kept_samples]

    st.subheader("Filtering Summary")
    a, b = st.columns([2, 1])

    with a:
        labels_values = [
            ("Selected submitters", len(selected_submitters)),
            ("Selected samples", len(selected_samples)),
            ("Samples after coverage filter", len(kept_samples)),
            ("Final rows for analysis", len(final_rows)),
        ]

        cols = st.columns(2)  # 2 metrics per row
        for i, (label, value) in enumerate(labels_values):
            with cols[i % 2]:
                st.metric(label, value)

    with b:
        if kept_samples:
            st.download_button(
                "Download kept samples",
                data="\n".join(sorted(kept_samples)),
                file_name="kept_samples.txt",
                mime="text/plain",
            )


    if len(final_rows) < 2:
        st.error("Insufficient data for PCA analysis. Try relaxing the filters.")
        return

    # PCA build + clean
    try:
        M, row_keys, col_keys = build_pca_matrix(final_rows, use_iupac, include_terminal)
        if M.shape[0] < 2:
            st.error("Need at least 2 comparisons for PCA")
            return

        M, col_keys = impute_and_clean_matrix(M, col_keys)
        if M.shape[1] == 0:
            st.error("No valid columns for PCA after preprocessing")
            return

        fig, pca_model = create_pca_plot(M, row_keys, show_spokes, show_centroids, point_size)
        st.subheader("PCA Visualization")
        st.plotly_chart(fig, use_container_width=True)

        st.subheader("Sample Contributions to Variance")
        var_df = calculate_variance_contributions(pca_model, col_keys)
        st.dataframe(var_df, use_container_width=True)
        st.download_button(
            "Download variance contributions",
            data=var_df.to_csv(index=False).encode("utf-8"),
            file_name="pca_variance_contributions.csv",
            mime="text/csv"
        )

        # Final filtered stats (full detail preserved)
        st.subheader("Final Filtered Statistics")
        f_sub = calculate_submitter_stats(final_rows, use_iupac, include_terminal)
        f_sam = calculate_sample_stats(final_rows, use_iupac, include_terminal)

        c1, c2 = st.columns(2)
        with c1:
            st.write("**Filtered Submitter Statistics**")
            disp = [
                "submitter", "samples", "rows", "valid", "pct_match",
                "pct_ins", "pct_ins_internal", "pct_ins_terminal",
                "pct_del", "pct_del_internal", "pct_del_terminal",
                "pct_sub", "ins_internal", "ins_terminal",
                "del_internal", "del_terminal"
            ]
            st.dataframe(f_sub[disp] if not f_sub.empty else pd.DataFrame(columns=disp),
                         use_container_width=True, hide_index=True)
            st.download_button(
                "Download filtered submitter stats",
                data=f_sub.to_csv(index=False).encode("utf-8"),
                file_name="filtered_submitter_stats.csv",
                mime="text/csv"
            )
        with c2:
            st.write("**Filtered Sample Statistics**")
            disp = ["sample", "rows", "valid", "pct_match", "pct_gap", 
                    "pct_gap_internal", "pct_gap_terminal", "pct_sub"]
            st.dataframe(f_sam[disp] if not f_sam.empty else pd.DataFrame(columns=disp),
                         use_container_width=True, hide_index=True)
            st.download_button(
                "Download filtered sample stats",
                data=f_sam.to_csv(index=False).encode("utf-8"),
                file_name="filtered_sample_stats.csv",
                mime="text/csv"
            )

    except Exception as e:
        st.error(f"Error in PCA analysis: {e}")
        st.error("Please check your data format and try again")


if __name__ == "__main__":
    main()