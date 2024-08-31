from mudata_explorer.app.url import load_url
import streamlit as st
from typing import Optional


def get_edit_view_flag() -> Optional[int]:
    val = st.query_params.get("edit-view")
    if val is not None:
        return int(val)


def set_edit_view_flag(ix: int):
    st.query_params["edit-view"] = ix


def clear_edit_view_flag():
    if "edit-view" in st.query_params:
        del st.query_params["edit-view"]


def get_editable_flag() -> bool:
    return st.query_params.get("editable", False)


def set_editable_flag(val: bool):
    assert isinstance(val, bool)
    if val:
        st.query_params["editable"] = val
    else:
        if "editable" in st.query_params:
            del st.query_params["editable"]


def update_edit_views():
    # Check the state of the input element used to update the editable flag
    flag = st.session_state["sidebar_toggle_editable"]
    # Make sure that the input element is in sync
    if flag != get_editable_flag():
        set_editable_flag(flag)


def check_file_url():
    if st.query_params.get("file"):
        load_url(st.query_params["file"])
        del st.query_params["file"]
