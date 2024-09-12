import streamlit as st
from mudata_explorer.helpers import save_load
from mudata_explorer.app.sidebar import setup_sidebar
from mudata_explorer.app.mdata import get_mdata_exists
from mudata_explorer.app.tables import build_mdata


def run():
    setup_sidebar()

    with st.container(border=1):
        st.write("#### Cirro Data Platform")

        st.page_link(
            "pages/cirro_load.py",
            label="Load Data From Cirro",
            icon=":material/download:"
        )

    save_load.upload_mdata(
        st.container(border=True)
    )

    # Let the user build a dataset by uploading tables
    with st.container(border=1):
        st.write("#### Upload Tables")
        st.page_link(
            "pages/load_tables.py",
            label="Build a MuData from Spreadsheets",
            icon=":material/upload:"
        )

    # Let the user browse existing datasets
    with st.container(border=1):
        st.write("#### Browse Public Datasets")
        st.page_link(
            "pages/public_data.py",
            label="Open Data from a Public Repository",
            icon=":material/public:"
        )
