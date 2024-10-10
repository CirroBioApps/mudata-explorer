from mudata_explorer.app.sidebar import setup_sidebar
from mudata_explorer.apps.microbiome import load
import streamlit as st

# Show the menu at the top
setup_sidebar(
    active_page="microbiome",
    title="Microbiome Report"
)

# Write the header
st.markdown('### Microbiome Report')

st.page_link(
    "https://living-figures.com/post/microbiome-report-intro/",
    label="What is a Microbiome Report?",
    icon=":material/info:",
    disabled=False
)

st.page_link(
    "pages/load_microbiome.py",
    label="Build a Microbiome Report",
    icon=":material/backup:",
    disabled=True
)

st.page_link(
    "https://living-figures.com/post/microbiome-report-gallery/",
    label="Browse the Gallery",
    icon=":material/photo_library:",
    disabled=False
)

load.run()
