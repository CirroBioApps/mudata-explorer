from mudata_explorer import app
import streamlit as st
from mudata_explorer.helpers import save_load, cirro


def run():
    app.setup_sidebar()

    # If the Cirro client has not been set up,
    # prompt the user to set it up
    if not st.session_state.get("Cirro"):
        cirro.setup_cirro_client(
            st.container(border=True)
        )

    # If it has been set up, show the save and load options
    else:

        cirro.load_from_cirro(
            st.container(border=True)
        )

        if app.has_mdata():
            cirro.save_to_cirro(
                st.container(border=True)
            )

    if app.has_mdata():
        save_load.download_mdata(
            st.container(border=True)
        )

    save_load.upload_mdata(
        st.container(border=True)
    )
