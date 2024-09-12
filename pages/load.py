import streamlit as st
from mudata_explorer.helpers import save_load
from mudata_explorer.app.sidebar import setup_sidebar


def run():
    setup_sidebar("load")

    st.subheader("Load MuData")

    with st.container(border=1):
        st.write("#### Cirro Data Platform")

        st.page_link(
            "pages/load_cirro.py",
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
            "pages/load_public_data.py",
            label="Open Data from a Public Repository",
            icon=":material/public:"
        )


if __name__ == "__main__":
    run()
