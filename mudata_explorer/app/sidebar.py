import streamlit as st
from streamlit.errors import StreamlitAPIException
from typing import Union, List
from mudata_explorer.app.mdata import get_mdata
from mudata_explorer.app.mdata import has_history
from mudata_explorer.app.query_params import get_editable_flag, update_edit_views
from mudata_explorer.app.query_params import check_file_url
from mudata_explorer.helpers.save_load import load_history
from mudata_explorer.helpers.plotting import plot_mdata
from streamlit.delta_generator import DeltaGenerator



def sidebar_page_links(page_links):
    for path, label, icon in page_links:
        st.sidebar.page_link(
            f"pages/{path}.py",
            label=label,
            icon=icon
        )


def sidebar_toggle_editable():

    st.sidebar.checkbox(
        "Edit Figures",
        value=get_editable_flag(),
        help="Display a set of menus to modify the figures.",
        on_change=update_edit_views,
        key="sidebar_toggle_editable"
    )


def sidebar_load_history():
    if get_mdata() is None:
        return

    if not has_history(exclude=['add_data', 'add_view']):
        return

    if st.sidebar.button("Rerun Analysis"):
        load_history()


def setup_sidebar(
    edit_views=False,
    load_history=False,
    page_layout="centered"
):
    """
    Setup the sidebar with links to all of the pages.
    If edit_views is True, add a checkbox to allow the user to edit the views.
    """
    try:
        st.set_page_config("MuData Explorer", layout=page_layout)
    except StreamlitAPIException:
        st.rerun()

    # If a file link is in the query params
    check_file_url()

    sidebar_page_links([
        ("save_load", "Save / Load", ":material/save:"),
        ("tables", "Tables", ":material/table:"),
        ("processes", "Analysis", ":material/function:"),
        ("views", "Figures", ":material/insert_chart:"),
        ("history", "History", ":material/history:"),
        ("about", "About", ":material/info:")
    ])
    if edit_views:
        sidebar_toggle_editable()
    if load_history:
        sidebar_load_history()

    plot_mdata(st.sidebar)


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
