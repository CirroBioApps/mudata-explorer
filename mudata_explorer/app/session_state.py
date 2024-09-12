import streamlit as st


def get_show_sidebar_flag() -> bool:
    return st.query_params.get("show_sidebar", "False") == 'True'


def set_show_sidebar_flag(val: bool):
    st.query_params["show_sidebar"] = str(val)
