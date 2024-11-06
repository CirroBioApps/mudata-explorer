from mudata_explorer.helpers import cirro
import streamlit as st


def run():
    # Let people load data from Cirro
    with st.container(border=1):
        st.write("#### Load From Cirro")
        cirro.load_from_cirro(
            filter_process_ids=[
                "process-hutch-virscan-1_0",
                "process-hutch-virscan-1_1",
                "process-hutch-virscan-1_2"
            ],
            show_link=False
        )
