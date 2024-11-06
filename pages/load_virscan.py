from mudata_explorer.app.sidebar import setup_sidebar
from mudata_explorer.apps.virscan import load
import streamlit as st

# Show the menu at the top
setup_sidebar(
    active_page="virscan",
    title="VirScan Report"
)

# Write the header
st.markdown('### VirScan Report')

load.run()
