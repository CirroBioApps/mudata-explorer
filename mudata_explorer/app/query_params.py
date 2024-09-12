from mudata_explorer.app.url import load_url
import streamlit as st
from typing import Optional


def get_edit_view_flag() -> Optional[int]:
    if "edit-view" in st.session_state:
        val = st.session_state["edit-view"]
        st.query_params["edit-view"] = val
        del st.session_state["edit-view"]
    else:
        val = st.query_params.get("edit-view")

    if val is not None:
        return int(val)


def set_edit_view_flag(ix: int):
    st.session_state["edit-view"] = ix
    st.switch_page("pages/view_details.py")


def clear_edit_view_flag():
    if "edit-view" in st.query_params:
        del st.query_params["edit-view"]


def check_file_url():
    if st.session_state.get("file"):
        url = st.session_state["file"]
        del st.session_state["file"]
    elif st.query_params.get("file"):
        url = st.query_params["file"]
        del st.query_params["file"]
    else:
        return

    load_url(url)
    st.switch_page("pages/views.py")
