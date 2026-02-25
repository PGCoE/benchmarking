import streamlit as st
from pathlib import Path
from typing import List, Tuple
import numpy as np
import pandas as pd
import plotly.express as px

import ui, metrics

#-------------------#
# Configuration
#-------------------# 

st.set_page_config(layout="wide")

#-------------------#
# Session State Initialization
#-------------------#

st.session_state.setdefault("warning_message", None)
st.session_state.setdefault("df", pd.DataFrame())
st.session_state.setdefault("df_final", pd.DataFrame())

st.session_state.progress = st.empty()
st.session_state.msg      = st.empty()
st.session_state.submsg   = st.empty()

st.session_state.fig_height = 500
st.session_state.fig_width  = 500

#-------------------#
# Sample Counter
#-------------------#
ui.sample_counter()

#-------------------#
# Side Bar
#-------------------#
with st.sidebar:
    st.header("PGCoE Viral Benchmarking Dashboard")

    with st.container(border=True):
        ui.load_data()

    if 'df_samplesheet' in st.session_state:
        ui.calc_metrics()

    if not st.session_state.df.empty:
        with st.container(border=True):
            if not st.session_state.df.empty:
                ui.data_processing_options()

        with st.container(border=True):
            if not st.session_state.df.empty:
                ui.apply_filters()

#-------------------#
# Data & Plots
#-------------------#
if not st.session_state.df.empty:

    tab1, tab2, tab3 = st.tabs(["Full Dataset", "Global View", "Sample View"])
        
    with tab1:
        st.dataframe(st.session_state.df)  
        st.info("💡 This table can be downloaded and re-used in place of the samplesheet!")        

    with tab2:
        if not st.session_state.df_final.empty:
            ui.global_view()
    
    with tab3:
        if not st.session_state.df_final.empty:
           ui.sample_view()  
        

elif st.session_state.warning_message:
    st.warning(st.session_state.warning_message)
else:
    st.session_state.msg.info("Upload data to get started.") 