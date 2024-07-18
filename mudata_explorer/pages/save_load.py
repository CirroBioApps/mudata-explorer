from mudata_explorer import app
import streamlit as st
from mudata_explorer.helpers import save_load


def run():
    app.setup_sidebar()

    with st.container(border=1):
        st.write("### Cirro: Save/Load")

        st.page_link(
            "pages/cirro_load.py",
            label="Load Data From Cirro",
            icon=":material/download:"
        )

        if app.has_mdata():
            st.page_link(
                "pages/cirro_save.py",
                label="Save Data To Cirro",
                icon=":material/save:"
            )

    if app.has_mdata():
        save_load.download_mdata(
            st.container(border=True)
        )

    save_load.upload_mdata(
        st.container(border=True)
    )
