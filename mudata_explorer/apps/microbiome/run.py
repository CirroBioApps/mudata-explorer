import streamlit as st
from streamlit.errors import StreamlitAPIException
from mudata_explorer.apps.microbiome import load
from mudata_explorer.apps.helpers import readme


def header():
    try:
        st.set_page_config(
            page_title='Microbiome Explorer',
            page_icon=':microbe:'
        )
    except StreamlitAPIException:
        st.rerun()
    st.title('Microbiome Explorer')
    st.write('Load microbiome data, generate plots, and explore the data.')


def run():

    header()
    load.run()
    st.write("---")
    st.markdown(readme("microbiome"))
