import streamlit as st


def get_show_sidebar_flag() -> bool:
    return st.session_state.get("show_sidebar", False)


def set_show_sidebar_flag(val: bool):
    st.session_state["show_sidebar"] = val
