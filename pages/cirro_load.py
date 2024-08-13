from mudata_explorer import app
import streamlit as st
from mudata_explorer.helpers import cirro


app.setup_sidebar()

st.write("#### Load from Cirro")
cirro.load_from_cirro()

st.page_link(
    'pages/save_load.py',
    label='Back',
    icon=":material/arrow_back:"
)
