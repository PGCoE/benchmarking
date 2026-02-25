import pandas as pd
import streamlit as st
import data_format 

FILT_COLS_EXCLUDE = set(['assembly', 'workflow'] + data_format.METRICS_COLS)

def get_counts(df, col_1, col_2, col):
    if col_1 not in df.columns or col_2 not in df.columns:
        return None

    df_subset = df[['sample', col_1, col_2]]

    count = (
        pd.concat([
            df_subset[['sample', col_1]],
            df_subset[['sample', col_2]].rename(columns={col_2: col_1}),
        ], ignore_index=True)
        .drop_duplicates()
        .groupby(col_1)
        .size()
        .reset_index(name='n')
        .sort_values(['n', col_1], ascending=[False, True])
        .reset_index(drop=True)
    )

    count[col] = count[col_1].astype(str) + " (" + count["n"].astype(str) + ")"
    return count


def apply_selection(df, col):
    col_1 = col + "_s1"
    col_2 = col + "_s2"

    count = get_counts(df, col_1, col_2, col)
    if count is None or count.empty:
        return {(col_1, col_2): []}
    
    opt_key = col + '_opts'
    pill_key      = col + "_pill"
    pill_key_prev = pill_key + "_prev"

    opts = ["All"] + count[col].tolist()
    if opt_key not in st.session_state:
        st.session_state[pill_key]      = ["All"]
        st.session_state[pill_key_prev] = ["All"]
    elif st.session_state[opt_key] != opts:
        st.session_state[pill_key]      = ["All"]
        st.session_state[pill_key_prev] = ["All"]

    st.session_state[opt_key] = opts
    
    if pill_key not in st.session_state:
        st.session_state[pill_key] = ["All"]
    if pill_key_prev not in st.session_state:
        st.session_state[pill_key_prev] = ["All"]

    def normalize():
        prev = st.session_state[pill_key_prev]
        cur  = st.session_state[pill_key]

        # newly added selections
        added = set(cur) - set(prev)

        if "All" in added:
            # user clicked All -> make it exclusive
            st.session_state[pill_key] = ["All"]
        elif "All" in cur and len(cur) > 1:
            # user added something else while All selected -> drop All
            st.session_state[pill_key] = [v for v in cur if v != "All"]
        elif not cur:
            # never allow empty -> default back to All
            st.session_state[pill_key] = ["All"]

        # IMPORTANT: update prev tracker after normalization
        st.session_state[pill_key_prev] = st.session_state[pill_key]
    
    with st.expander(col):
        st.pills(None, opts, selection_mode="multi", key=pill_key, on_change=normalize)

    select = st.session_state[pill_key]
    select_orig = "All" if "All" in select else count.loc[count[col].isin(select), col_1].tolist()
    return {(col_1, col_2): select_orig}

