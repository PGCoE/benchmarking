import numpy as np
import pandas as pd
import plotly.express as px
import streamlit as st

from scipy.spatial.distance import squareform
import scipy.cluster.hierarchy as sch
import plotly.express as px
import plotly.graph_objects as go
import plotly.figure_factory as ff
from plotly.subplots import make_subplots

from Bio.Align import PairwiseAligner

import io_ops
import metrics

def classical_pcoa(D, k=2):
    """
    Classical PCoA from a (n x n) distance matrix D.
    Returns:
      coords: (n x k) coordinates
      var: (k,) fraction variance explained
      eigvals: all eigenvalues (for diagnostics)
    """
    D = np.asarray(D, dtype=float)
    n = D.shape[0]

    # Double-center the squared distances
    J = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * J @ (D ** 2) @ J

    # Eigen decomposition (symmetric)
    eigvals, eigvecs = np.linalg.eigh(B)

    # Sort descending
    idx = np.argsort(eigvals)[::-1]
    eigvals = eigvals[idx]
    eigvecs = eigvecs[:, idx]

    # Keep only positive eigenvalues (numerical tolerance)
    tol = 1e-12
    pos = eigvals > tol
    eigvals_pos = eigvals[pos]
    eigvecs_pos = eigvecs[:, pos]

    # Compute coords
    k = min(k, eigvals_pos.size)
    coords = eigvecs_pos[:, :k] * np.sqrt(eigvals_pos[:k])

    # Variance explained
    var = eigvals_pos[:k] / eigvals_pos.sum()

    return coords, var, eigvals

def plot_pcoa():
    D = st.session_state.matrix_values
    labels = st.session_state.matrix_labels

    coords, var, eigvals = classical_pcoa(D, k=2)

    try:
        df_plot = pd.DataFrame(coords, columns=["PCo1", "PCo2"])
        df_plot["label"] = labels
    except:
        st.session_state.fig_pcoa = None
        return

    xlab = f"PCo1 ({var[0]*100:.1f}%)"
    ylab = f"PCo2 ({var[1]*100:.1f}%)" if len(var) > 1 else "PCo2"

    fig = px.scatter(
        df_plot,
        x="PCo1",
        y="PCo2",
        text="label",
        hover_name="label",
        hover_data=df_plot.columns.tolist(),
        title="PCoA (classical) on Pairwise Dissimilarities",
    )

    fig.update_traces(textposition="top center", marker=dict(size=10))
    fig.update_layout(
        xaxis_title=xlab,
        yaxis_title=ylab,
        width=st.session_state.fig_width,
        height=st.session_state.fig_height,
        margin=dict(l=40, r=40, t=60, b=60),
    )

    st.session_state.fig_pcoa = fig


def plot_heatmap(method="average", title="Pairwise Dissimilarity Heatmap"):
    # Pull from session state — same source as plot_pcoa
    D      = st.session_state.matrix_values
    labels = st.session_state.matrix_labels

    mat = pd.DataFrame(D, index=labels, columns=labels)

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
    fig.update_yaxes(autorange="reversed", row=1, col=1)
    fig.update_xaxes(autorange="reversed", row=1, col=1)
    fig.update_xaxes(showticklabels=False, row=1, col=1)
    fig.update_yaxes(visible=False, row=1, col=1)

    fig.update_xaxes(side="bottom", row=2, col=1)
    fig.update_yaxes(autorange="reversed", row=2, col=1)

    fig.update_layout(
        width=st.session_state.fig_width,
        height=st.session_state.fig_height,
        title=title,
        margin=dict(l=10, r=10, t=40, b=10),
    )

    st.session_state.fig_heatmap = fig

def plot_msa(sample, s1, s2):

    if not sample or not s1 or not s2:
        return
    
    df = st.session_state.df_final

    df_subset = df[
        (df["sample"] == sample) &
        (df["source_s1"] == s1) &
        (df["source_s2"] == s2)
    ]
    # try other way - sometimes group variables are swapped
    if df_subset.empty:
        df_subset = df[
        (df["sample"] == sample) &
        (df["source_s1"] == s2) &
        (df["source_s2"] == s1)
    ]
        
    if df_subset.empty:
        st.error("No rows matched sample/source filters.")
        return

    row = df_subset.iloc[0]

    try:
        seq1 = io_ops.load_one_sequence(row["assembly_s1"])
        seq2 = io_ops.load_one_sequence(row["assembly_s2"])
    except Exception as e:
        st.exception(e)
    aln, revcomp = metrics.align_pair(seq1, seq2)

    match, valid, invalid, ins, del_, sub = df_subset[["match", "valid", "invalid", "ins", "del", "sub"]].iloc[0]
    pid = 100 * match / valid
    overlap = 100 * valid / (valid + invalid)
    termini_status = 'invalid' if st.session_state.ignore_termini else 'valid'
    target = s1 + " (reverse complement)" if revcomp else s1
    query = s2

    header = [
        f"Sample: {sample}",
        f"Target: {target}",
        f"Query: {query}",
        "",
        f"Overlap: {overlap:.2f}%",
        f"Identity: {pid:.2f}%",
        "",
        f"valid: {valid}",
        f"invalid: {invalid}",
        f"ins {ins}",
        f"del: {del_}",
        f"sub: {sub}",
        f"match: {match}",
        "",
        f"** Invalid sites not counted in identity calculation **",
        f"** Terminal indels treated as {termini_status} sites **",
        "",
        ""
    ]

    alignment_text = '\n'.join(header) + str(aln)

    st.code(alignment_text)

def plot_dist(group=False):
    df_final = st.session_state.df_final
    if df_final.empty:
        return

    gcol = st.session_state.get('gcol', '')
    gcol_1, gcol_2 = gcol + "_s1", gcol + "_s2"

    if group:
        df_swapped = df_final.copy()
        df_swapped[gcol_1], df_swapped[gcol_2] = df_final[gcol_2].values, df_final[gcol_1].values
        df = pd.concat([df_final, df_swapped], ignore_index=True)
    else:
        df = df_final

    custom_data = ["sample", gcol_1, gcol_2, "percent_overlap"]
    if gcol != "source":
        custom_data += ["source_s1", "source_s2"]

    x_col = gcol_1 if group else "sample"
    category_order = (
        df.groupby(x_col)["percent_match"]
        .mean()
        .sort_values()
        .index.tolist()
    )

    fig = px.box(
        df,
        x=x_col,
        y="percent_match",
        points="all",
        custom_data=custom_data,
        category_orders={x_col: category_order}
    )

    fig.update_traces(hovertemplate=(
        "sample: %{customdata[0]}<br>"
        f"{gcol_1}: %{{customdata[1]}}<br>"
        f"{gcol_2}: %{{customdata[2]}}<br>"
        "match: %{y}%<br>"
        f"overlap: %{{customdata[3]:.2f}}%<br>"
        "<extra></extra>"
    ))

    return fig