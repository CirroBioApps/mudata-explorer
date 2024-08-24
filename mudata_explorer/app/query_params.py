import streamlit as st
from mudata_explorer.app.url import load_url


def get_edit_views_flag() -> bool:
    return st.query_params.get("editable", False)


def set_edit_views_flag(val: bool):
    assert isinstance(val, bool)
    if val:
        st.query_params["editable"] = val
    else:
        if "editable" in st.query_params:
            del st.query_params["editable"]


def update_edit_views():
    flag = st.session_state["sidebar_edit_views"]
    if flag != get_edit_views_flag():
        set_edit_views_flag(flag)


def check_file_url():
    if st.query_params.get("file"):
        load_url(st.query_params["file"])
        del st.query_params["file"]
