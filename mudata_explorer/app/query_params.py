from copy import copy
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
    if st.query_params.get("file"):
        load_url(st.query_params["file"])
        del st.query_params["file"]
        st.switch_page("pages/views.py")
