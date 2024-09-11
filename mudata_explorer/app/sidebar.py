import streamlit as st
from streamlit.errors import StreamlitAPIException
from typing import Union, List
from mudata_explorer.app import mdata
from mudata_explorer.app.query_params import get_show_sidebar_flag, set_show_sidebar_flag
from mudata_explorer.app.query_params import get_show_menu_flag
from mudata_explorer.app.query_params import check_file_url
from mudata_explorer.helpers.save_load import load_history
from streamlit.delta_generator import DeltaGenerator



def sidebar_page_links(page_links):
    for path, label, icon in page_links:
        st.sidebar.page_link(
            f"pages/{path}.py",
            label=label,
            icon=icon
        )


def sidebar_toggle_button():

    is_open = get_show_sidebar_flag()

    if st.sidebar.button(
        f"{'Close' if is_open else 'Open'} Sidebar",
        key="sidebar_toggle_button"
    ):
        st.query_params['show_sidebar'] = str(not is_open)
        st.session_state['show_sidebar'] = str(not is_open)
        st.rerun()


def sidebar_load_history(id="main"):
    if not mdata.get_mdata_exists(id=id):
        return

    if not mdata.has_history(exclude=['add_data', 'add_view'], id=id):
        return

    if st.sidebar.button("Rerun Analysis"):
        load_history()


def setup_sidebar(
    sidebar_toggle=False,
    load_history=False,
    page_layout="centered"
):
    """
    Setup the sidebar with links to all of the pages.
    If sidebar_toggle is True, add a checkbox to allow the user to edit the views.
    """
    try:
        st.set_page_config(
            "MuData Explorer",
            layout=page_layout,
            initial_sidebar_state="expanded" if get_show_menu_flag() else "collapsed"
        )
    except StreamlitAPIException:
        st.rerun()

    # If a file link is in the query params
    check_file_url()

    # Check if there are elements in the session state
    # that need to be propagated to the query params
    if 'show_sidebar' in st.session_state:
        set_show_sidebar_flag(st.session_state['show_sidebar'])

    sidebar_page_links([
        ("save_load", "Save / Load", ":material/save:"),
        ("tables", "Tables", ":material/table:"),
        ("processes", "Analysis", ":material/function:"),
        ("views", "Figures", ":material/insert_chart:"),
        ("history", "History", ":material/history:"),
        ("public_data", "Public Data", ":material/local_library:"),
        ("about", "About", ":material/info:")
    ])
    if sidebar_toggle:
        sidebar_toggle_button()
    if load_history:
        sidebar_load_history()


def landing_shortcuts():

    show_shortcuts([
        ("tables", "Upload Tables (*.csv)", ":material/table:"),
        ("about", "About", ":material/info:")
    ])


def show_shortcuts(
    shortcuts: List[tuple],
    ncol=3,
    container: Union[None, DeltaGenerator] = None
):
    cols = None

    for ix, (path, label, icon) in enumerate(shortcuts):
        if ix % ncol == 0:
            cols = (
                st.columns(ncol)
                if container is None
                else container.columns(ncol)
            )

        cols[ix % 3].page_link(
            f"pages/{path}.py",
            label=label,
            icon=icon,
            use_container_width=True
        )
