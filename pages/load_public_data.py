from mudata_explorer.app.sidebar import setup_sidebar
from mudata_explorer.public_data import repositories
import streamlit as st


def run():

    setup_sidebar("load_public_data")

    st.subheader("Load MuData: From Public Data")

    st.markdown(
        """
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

    st.page_link(
        "pages/load.py",
        label="Back",
        icon=":material/arrow_back:"
    )


if __name__ == "__main__":
    run()
