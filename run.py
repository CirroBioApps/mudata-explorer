#!/usr/bin/env streamlit run

from mudata_explorer import app
import streamlit as st

app.setup_pages()
st.write("### MuData Explorer")
app.landing_shortcuts()
