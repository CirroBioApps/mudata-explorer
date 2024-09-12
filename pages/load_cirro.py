from mudata_explorer.app.sidebar import setup_sidebar
import streamlit as st
from mudata_explorer.helpers import cirro


setup_sidebar("load_cirro")

st.subheader("Load MuData: From Cirro")

cirro.load_from_cirro()

st.page_link(
    'pages/load.py',
    label='Back',
    icon=":material/arrow_back:"
)
