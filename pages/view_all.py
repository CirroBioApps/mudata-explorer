from mudata_explorer.app.sidebar import setup_sidebar
from mudata_explorer.app.mdata import get_mdata_exists
import streamlit as st
from mudata_explorer.app.mdata import get_views, get_mdata_exists
from mudata_explorer.helpers.assets import make_view
from mudata_explorer.app.session_state import set_show_sidebar_flag


def make_views():

    views = get_views()

    return [
        make_view(
            ix=ix,
            type=view["type"],
            params=view["params"]
        )
        for ix, view in enumerate(views)
    ]


def view_non_editable():

    # All of the views defined in the dataset
    mdata_views = make_views()

    for view in mdata_views:

        # Attach the view to the display
        view.attach()


def run():

    if not get_mdata_exists():
        st.switch_page("pages/load.py")

    set_show_sidebar_flag(False)

    setup_sidebar("view_all")

    # Show the views with no edit buttons
    view_non_editable()


if __name__ == "__main__":
    run()
