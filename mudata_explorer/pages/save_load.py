from mudata_explorer import app
import streamlit as st
from streamlit.delta_generator import DeltaGenerator
from mudata_explorer.helpers import save_load


def run():
    app.setup_sidebar()

    if app.has_mdata():
        save_load.download_mdata(
            st.container(border=True)
        )

    save_load.upload_mdata(
        st.container(border=True)
    )
