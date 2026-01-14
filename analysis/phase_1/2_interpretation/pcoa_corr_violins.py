#! /usr/bin/env python3

"""
Tool for parsing benchmark viral consensus assembly results from TheiaViral/TheiaCoV (Public Health Bioinformatics task_consensus_qc) and VAPER
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from collections import defaultdict

def compile_corr_csvs(csv_paths):
    """Read correlation CSV to file""" 
    df = pd.DataFrame()
    for csv_path in csv_paths:
        try:
            t_df = pd.read_csv(csv_path, index_col=False)
        except pd.errors.EmptyDataError:
            print(f"Warning: {csv_path} is empty, skipping...")
            continue
        t_df = t_df.replace(["NA"], 0)
        t_df = t_df.astype({"pco1_r2": float, "pco2_r2": float})
        t_df["tot"] = t_df["pco1_r2"] + t_df["pco2_r2"]
        df = pd.concat([df, t_df[["var", "tot"]]])
    return df


def make_plots(tot_df, title, output="pcoa_axes_correlation.html"):
    """Make violin plots summarizing benchmarking results"""

    # initialize figure
    fig = go.Figure()

    # add violin
    fig.add_trace(go.Violin(x=tot_df["var"],
                            y=tot_df["tot"]))
#                            legendgroup=, scalegroup=software, name=software,
 #                           line_color=colors[software]))

    fig.update_yaxes(range=[0, 1])
    fig.update_traces(box_visible=True, meanline_visible=True, points="outliers")
    fig.update_layout(violinmode="group",
                      yaxis_title="Total R² with PCoA Axes 1 & 2",
                      title_text=title,
                      font=dict(size=20, color="black"))
    fig.write_html(output)



def main(csv_paths, title, output):

    tot_df = compile_corr_csvs(csv_paths)
    if not output:
        output = "pcoa_axes_correlation.html"

    make_plots(tot_df, title, output) 


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse PGCoE PCoA correlation CSVs and make violin plot")
    parser.add_argument("-c", "--csv", nargs='*', help="CSV paths")
    parser.add_argument("-t", "--title", help="Chart title")
    parser.add_argument("-o", "--output", help="Output HTML name")
    args = parser.parse_args()

    main(args.csv, args.title, args.output)
    sys.exit(0)