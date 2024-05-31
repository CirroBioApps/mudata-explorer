from mudata_explorer import app
import streamlit as st


def run():
    app.setup_sidebar()

    # Read the contents of the Readme.md file
    with open("README.md", "rt") as handle:
        readme = handle.read()

    st.write(readme)
