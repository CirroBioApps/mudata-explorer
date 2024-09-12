import streamlit as st
from mudata_explorer.app.sidebar import setup_sidebar


def run():
    setup_sidebar("about")

    # Read the contents of the Readme.md file
    with open("README.md", "rt") as handle:
        readme = handle.read()

    st.write(readme)


if __name__ == "__main__":
    run()
