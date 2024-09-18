from mudata_explorer.app.sidebar import setup_sidebar
from mudata_explorer.app.mdata import get_mdata_exists
import streamlit as st
from mudata_explorer.app.mdata import get_mdata_exists
from mudata_explorer.app.session_state import set_show_sidebar_flag
from mudata_explorer.app.view import view_non_editable


def run():

    if not get_mdata_exists():
        st.switch_page("pages/load.py")

    set_show_sidebar_flag(False)

    setup_sidebar("view_all")

    # Show the views with no edit buttons
    view_non_editable()


if __name__ == "__main__":
    run()
