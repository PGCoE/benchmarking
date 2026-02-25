import pandas as pd
import streamlit as st
from pathlib import Path

import io_ops
import data_format
import data_filter
import data_plot
import metrics

def sample_counter():
    st.session_state.setdefault("samples_f0", 0)
    st.session_state.setdefault("samples_f1", 0)
    st.session_state.setdefault("samples_f2", 0)
    st.session_state.setdefault("samples_fi", 0)


    with st.container():
        st.caption(f"Sample counts: Total={st.session_state['samples_f0']}, Filter_1={st.session_state['samples_f1']}, Filter_2={st.session_state['samples_f2']}, Final={st.session_state['samples_fi']}")
    
def load_data():
    st.header("Load Samplesheet")
    samplesheet_path = st.file_uploader(
        "Samplesheet",
        type=["csv"],
        accept_multiple_files=False
    )

    # Only reset and reprocess when files are actually provided
    if not samplesheet_path:
        # No files uploaded yet — leave session state untouched
        return
    
    df = io_ops.load_csv(samplesheet_path)
    if set(data_format.METRICS_COLS).issubset(df.columns):
        st.session_state.df = df
    elif set(io_ops.SAMPLESHEET_COLS).issubset(df.columns):
        st.session_state.df_samplesheet = df
    else:
        st.error(f"Input is missing required columns")

def data_processing_options():
    st.header("Data Processing Options")

    df = st.session_state.df

    # Grouping variable
    st.markdown("### Sample grouping variable")
    opts = ['source'] 
    for c in df.columns:
        if c in data_format.METRICS_COLS:
            continue
        if c.endswith('_s1') or c.endswith('_s2'):
            c = c[:-3]
        if c not in opts:
            opts.append(c)
        
    gcol = st.selectbox("Group by (default: `source`)", opts)
    if gcol:
        data_format.apply_grouping(gcol, 'gcol')

    st.session_state.gcol = gcol
    st.session_state.g_opts = opts

    # Terminal indels
    st.markdown("### Terminal indels")
    ignore_termini = st.toggle('Ignore Termini', value=False, help="Insertions / Deletions at the end of sequences will be considered `invalid` and not included in the similarity / dissimilarity calculations.")

    # Add row stats
    st.session_state.df_stats = data_format.row_stats(df, ignore_termini)


def apply_filters():
    df     = st.session_state.df_stats
    gcol   = st.session_state.gcol
    g_opts = st.session_state.g_opts

    st.session_state.df_final = pd.DataFrame()

    if df.empty:
        return
    
    st.session_state["samples_f0"] = len(df)

    st.header("Apply filters")
    
    st.markdown("### Filter 1: Minimum sequence overlap")
    min_valid = st.number_input(
        "Select a threshold (default: `0.80`)",
        min_value=0.0,
        max_value=1.0,
        step=0.01,
        value=0.80,
        format="%.2f",
        help=f"Minimum overlap (valid_frac) in sample assemblies between two {gcol}s to be considered"
    )
 
    df_subset = df[df["valid_frac"] > min_valid]
    st.session_state["samples_f1"] = len(df_subset)

    st.markdown("### Filter 2: Include variables")
    filt_cols  = [ col for col in set(g_opts) if col not in data_filter.FILT_COLS_EXCLUDE ]
    selections = []
    for col in filt_cols:
        selections.append(data_filter.apply_selection(df_subset, col))

    if selections:
        valid_cols = df_subset.columns
        for rec in selections:
            for (col_1, col_2), values in rec.items():
                if col_1 not in valid_cols or col_2 not in valid_cols:
                    continue
                if 'All' in values:
                    continue
                df_subset = df_subset[
                    df_subset[col_1].isin(values) & df_subset[col_2].isin(values)
                ]

    if df_subset.empty:
        st.session_state.msg.warning("No data remains after filtering. Try relaxing your filter thresholds.")
        return

    sample_map = {}
    for rec in df_subset[['sample', 'gcol_a', 'gcol_b']].to_dict(orient='records'):
        sample = rec['sample']
        sample_map.setdefault(rec['gcol_a'], set()).add(sample)
        sample_map.setdefault(rec['gcol_b'], set()).add(sample)

    sample_sets = list(sample_map.values())

    if not sample_sets:
        st.session_state.msg.warning("No grouping data found after filtering. Check your grouping variable selection.")
        return

    sample_inter = set.intersection(*sample_sets)

    if not sample_inter:
        st.session_state.msg.warning(f"No samples are common across all {gcol} groups after filtering.")
        return

    df_final = df_subset[df_subset['sample'].isin(sample_inter)]

    st.session_state["samples_f2"] = len(df_final)
    st.session_state["samples_fi"] = len(df_final)

    st.session_state.sample_inter = sample_inter
    st.session_state.df_final     = df_final


def global_view():
    st.session_state.fig_heatmap     = None
    st.session_state.fig_pcoa        = None

    st.caption(f"Filtered Samples: {len(st.session_state.sample_inter)}")

    data_format.global_stats()
    st.subheader("Global Statistics (After Filtering)")
    st.dataframe(st.session_state.df_global)


    data_format.to_distance_matrix()
    data_plot.plot_pcoa()
    data_plot.plot_heatmap()

    col1, col2 = st.columns(2)
    with col1:
        with st.container(border=True): 
            fig_heatmap = st.session_state.fig_heatmap
            if fig_heatmap:
                st.plotly_chart(fig_heatmap)
            else:
                st.warning("Heatmap could not be generated")

    with col2:
        with st.container(border=True): 
            fig_pcoa = st.session_state.fig_pcoa
            if fig_pcoa:
                st.plotly_chart(fig_pcoa)
            else:
                st.warning("PCoA could not be generated")

def sample_view():
    st.session_state.fig_sample_dist = None
    data_plot.plot_sample_dist()

    fig_sample_dist = st.session_state.fig_sample_dist
    if fig_sample_dist:
        sample_selection = st.plotly_chart(fig_sample_dist, use_container_width=False, config={"displayModeBar": False}, on_select='rerun', selection_mode='points')

        points = sample_selection.get('selection', {}).get('points', [])

        if not points:
            st.info("Select a data point to display the multiple sequence alignment")
        else:
            point0 = points[0]
            s1, s2 = point0.get('customdata', [None, None])
            sample = point0.get('x', None)

            if sample and s1 and s2:
                st.session_state.sample_selection = (sample, s1, s2)
                st.divider()
                data_plot.plot_msa()
    else:
        st.warning("Sample distributions could not be generated")

def calc_metrics():
    with st.container(border=True):
        st.header("Calculate Metrics")

        outdir = st.text_input("Output Directory", "results/")
        st.session_state.outdir = Path(outdir)

        if st.button("Calculate"):
            st.toast(
                "Metrics can be downloaded from the 'Full Dataset' tab and uploaded in place of the samplesheet",
                icon = "💡"
            ) 
            metrics.run()