from mudata_explorer.app.sidebar import setup_sidebar
import streamlit as st
from mudata_explorer.helpers import cirro


setup_sidebar()

st.write("#### Load from Cirro")
cirro.load_from_cirro()

st.page_link(
    'pages/save_load.py',
    label='Back',
    icon=":material/arrow_back:"
)
