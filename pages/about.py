from mudata_explorer import app
import streamlit as st

if __name__ == "__main__":
    app.setup_pages()

    # Read the contents of the Readme.md file
    with open("README.md", "rt") as handle:
        readme = handle.read()

    st.write(readme)
