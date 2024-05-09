#!/usr/bin/env streamlit run

from mudata_explorer import app
import streamlit as st

app.setup_pages()
st.write("""
### MuData Explorer

Welcome to the MuData Explorer! This app is designed to help you visualize
and interact with your data.

- **Upload Tables**: Add data to the app by uploading CSV files.
- **Run Analysis**: Apply statistical analyses and save the results.
- **View Data**: Explore your data with tables and plots.
- **Save / Load Dataset**: Download a snapshot of your dataset to share with others.

""")

app.landing_shortcuts()
