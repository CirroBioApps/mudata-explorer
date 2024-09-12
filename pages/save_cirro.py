from mudata_explorer.app.sidebar import setup_sidebar
from mudata_explorer.helpers import cirro, save_load
from mudata_explorer.app.mdata import get_mdata_exists
import streamlit as st

def run():
    if not get_mdata_exists():
        st.switch_page("pages/load.py")

    setup_sidebar("save_cirro")

    st.write("### Save MuData: Cirro Data Platform")

    cirro.save_to_cirro()

    st.page_link(
        'pages/save.py',
        label='Back',
        icon=":material/arrow_back:"
    )


if __name__ == "__main__":
    run()