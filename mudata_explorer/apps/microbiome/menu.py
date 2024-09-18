from mudata_explorer.app.sidebar import setup_sidebar
import streamlit as st


def show_menu(active_page: str):

    setup_sidebar(
        active_page="microbiome",
        title="Microbiome Report"
    )

    # Write the header
    st.markdown('### Microbiome Report')

    st.page_link(
        "pages/microbiome.py",
        label="What is a Microbiome Report?",
        icon=":material/info:",
        disabled=active_page == "microbiome"
    )

    st.page_link(
        "pages/load_microbiome.py",
        label="Build a Microbiome Report",
        icon=":material/backup:",
        disabled=active_page == "load_microbiome"
    )

    st.page_link(
        "pages/gallery_microbiome.py",
        label="Browse the Gallery",
        icon=":material/photo_library:",
        disabled=active_page == "gallery_microbiome"
    )
