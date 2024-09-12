import streamlit as st
from mudata_explorer.app.tables import build_mdata
from mudata_explorer.app.sidebar import setup_sidebar


def run():

    setup_sidebar("load_tables")

    st.subheader("Load MuData: From Tables")

    # Let the user build a dataset by uploading tables
    with st.container(border=1):
        st.write("### Upload Tables")
        build_mdata()

    st.page_link(
        'pages/load.py',
        label='Back',
        icon=":material/arrow_back:"
    )


if __name__ == "__main__":
    run()
