from mudata_explorer.app.sidebar import setup_sidebar
from mudata_explorer.helpers import cirro, save_load
from mudata_explorer.app.mdata import get_mdata_exists
import streamlit as st

def run():
    if not get_mdata_exists():
        st.switch_page("pages/load.py")

    setup_sidebar("save")

    st.page_link(
        'pages/save_cirro.py',
        label='Save to Cirro',
        icon=":material/public:"
    )

    save_load.download_mdata(
        st.container(border=True)
    )


if __name__ == "__main__":
    run()