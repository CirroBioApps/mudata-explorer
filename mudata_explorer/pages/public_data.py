from mudata_explorer.app.sidebar import setup_sidebar
from mudata_explorer.public_data import repositories
import streamlit as st


def run():

    setup_sidebar()

    st.markdown(
        """
        # Public Data

        This page contains a collection of public data repositories that can be used for analysis.

        """
    )

    # Let the user select which repository to view
    repo_name = st.selectbox(
        "Select a repository",
        list(repositories.keys())
    )

    # Run that method
    repositories[repo_name].run()
