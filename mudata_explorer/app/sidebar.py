import streamlit as st
from streamlit.errors import StreamlitAPIException
from typing import Union, List
from mudata_explorer.app import mdata
from mudata_explorer.app.query_params import check_file_url
from mudata_explorer.helpers.save_load import load_history
from mudata_explorer.app.mdata import get_mdata_exists
from streamlit.delta_generator import DeltaGenerator



def sidebar_page_links(page_links, disabled_pages=[]):
    for path, label, icon in page_links:
        st.sidebar.page_link(
            f"pages/{path}.py",
            label=label,
            icon=icon,
            disabled=path in disabled_pages
        )


def sidebar_load_history(id="main"):
    if not mdata.get_mdata_exists(id=id):
        return

    if not mdata.has_history(exclude=['add_data', 'add_view'], id=id):
        return

    if st.sidebar.button("Rerun Analysis"):
        load_history()


def setup_sidebar(active_page: str, title="MuData Explorer"):
    """
    Setup the sidebar with links to all of the pages.
    If sidebar_toggle is True, add a checkbox to allow the user to edit the views.
    """
    wide_pages = ["view_sidebar", "view_details", "processes"]
    hide_menu_pages = ["view_all", "view_sidebar", "microbiome"]
    try:
        st.set_page_config(
            title,
            layout="wide" if active_page in wide_pages else "centered",
            initial_sidebar_state="collapsed" if active_page in hide_menu_pages else "expanded"
        )
    except StreamlitAPIException:
        st.rerun()

    if active_page.startswith("view_"):
        # If we navigate back to the Figures page, use this one
        st.session_state["last-view-page"] = active_page

    # If a file link is in the query params
    check_file_url()

    # If there is no data available
    if not get_mdata_exists():
        # Disable the pages that show data
        disabled_pages = [
            "tables",
            "processes",
            "view_all",
            "history",
            "save"
        ]
    else:
        disabled_pages = []

    sidebar_page_links(
        [
            ("load", "Load", ":material/backup:"),
            ("save", "Save", ":material/save:"),
            ("tables", "Tables", ":material/table:"),
            ("processes", "Analysis", ":material/function:"),
            (st.session_state.get("last-view-page", "view_all"), "Figures", ":material/insert_chart:"),
            ("history", "History", ":material/history:"),
            ("about", "About", ":material/info:")
        ],
        disabled_pages=disabled_pages
    )
    if active_page in ["view_sidebar", "view_details", "view_all", "microbiome"]:
        st.sidebar.write('---')
        st.sidebar.markdown("**View Mode**")
        sidebar_page_links(
            [
                ("view_all", "Figures Only", ":material/insert_chart:"),
                ("view_sidebar", "Figures + Sidebar", ":material/edit:")
            ],
            disabled_pages=[active_page]
        )

    elif active_page == "processes":
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
