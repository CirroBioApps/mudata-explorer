from copy import copy
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


def get_show_menu_flag() -> bool:
    return st.query_params.get("show_menu", True)


def set_show_menu_flag(val: bool):
    assert isinstance(val, bool)
    if not val:
        st.query_params["show_menu"] = val
    else:
        if "show_menu" in st.query_params:
            del st.query_params["show_menu"]


def toggle_show_menu():
    set_show_menu_flag(not get_show_menu_flag())


def get_show_sidebar_flag() -> bool:
    return st.query_params.get("show_sidebar", "False") == 'True'


def set_show_sidebar_flag(val: bool):
    st.query_params["show_sidebar"] = str(val)


def toggle_show_sidebar():
    set_show_sidebar_flag(not get_show_sidebar_flag())


def check_file_url():
    if st.query_params.get("file"):
        load_url(st.query_params["file"])
        del st.query_params["file"]
        st.switch_page("pages/views.py")
