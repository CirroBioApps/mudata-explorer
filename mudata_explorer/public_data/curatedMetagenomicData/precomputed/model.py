from mudata_explorer.public_data.curatedMetagenomicData.inventory import long_description
import streamlit as st


def run():
    
    # Show the precomputed inventory
    st.markdown(long_description)