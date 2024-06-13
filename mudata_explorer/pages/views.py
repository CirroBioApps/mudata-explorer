from mudata_explorer import app
from mudata_explorer.base.view import View
from mudata_explorer.helpers import all_views, make_view
from mudata_explorer.helpers import asset_categories, asset_dataframe
from mudata_explorer.helpers import filter_by_category
import streamlit as st
from streamlit.delta_generator import DeltaGenerator


def make_views(editable=False):

    views = app.get_views()

    return [
        make_view(
            ix=ix,
            editable=editable,
            **view
        )
        for ix, view in enumerate(views)
    ]


def make_editable(ix: int):
    st.query_params["edit-view"] = ix


def edit_view(view: View, container: DeltaGenerator, ix: int, n_views: int):

    # If the views are not editable, don't show the edit buttons
    if not app.get_edit_views_flag():
        # Instead, just set up the params for the view
        view.get_data()
        return

    cols = container.columns([1, 1, 1, 1, 1, 5])

    # The first button allows the user to provide inputs
    cols[0].button(
        ":pencil:",
        help="Edit the inputs for this view.",
        key=f"edit-inputs-{ix}",
        on_click=make_editable,
        args=(ix,)
    )

    # For every view past the first
    cols[1].button(
        ":arrow_up_small:",
        on_click=move_up,
        key=f"move-up-{ix}",
        args=(ix,),
        help="Move this view up in the list.",
        disabled=ix == 0
    )

    if cols[2].button(
        ":heavy_minus_sign:",
        key=f"delete-view-{ix}",
        args=(ix,),
        help="Delete this view."
    ):
        app.delete_view(ix)
        st.rerun()
    cols[3].button(
        ":heavy_plus_sign:",
        on_click=app.duplicate_view,
        key=f"duplicate-view-{ix}",
        args=(ix,),
        help="Duplicate this view."
    )
    cols[4].button(
        ":arrow_down_small:",
        on_click=move_down,
        key=f"move-down-{ix}",
        args=(ix,),
        help="Move this view down in the list.",
        disabled=ix == (n_views - 1)
    )


def button_add_view():

    st.write("#### Add a new view")

    # Let the user select the type of view to add
    all_categories = asset_categories(all_views)

    selected_category = st.selectbox(
        "Select a category",
        all_categories
    )

    # Let the user select which view to add, filtering by category
    filtered_views = filter_by_category(all_views, selected_category)

    # Get the assets needed to select from the filtered views
    df = asset_dataframe(filtered_views)

    selected_name = st.selectbox(
        "Select a view to add",
        df["name"].tolist()
    )
    selected_ix = df["name"].tolist().index(selected_name)
    selected_type = df["type"].tolist()[selected_ix]
    help_text = df["help_text"].tolist()[selected_ix]

    if help_text is not None:
        st.markdown(help_text)

    # Instantiate the selected view type if the user clicks a button
    st.button(
        f"Add {selected_name}",
        on_click=button_add_view_callback,
        args=(selected_type,),
        use_container_width=True
    )


def button_add_view_callback(selected_type: str):
    app.add_view(selected_type)
    st.query_params["edit-view"] = len(app.get_views()) - 1


def move_up(ix: int):
    views = app.get_views()
    (
        views[ix - 1],
        views[ix]
    ) = (
        views[ix],
        views[ix - 1]
    )
    app.set_views(views)


def move_down(ix: int):
    views = app.get_views()
    (
        views[ix + 1],
        views[ix]
    ) = (
        views[ix],
        views[ix + 1]
    )
    app.set_views(views)


def run():

    app.setup_sidebar(
        edit_views=True,
        page_layout=(
            "centered"
            if st.query_params.get("edit-view") is None
            else "wide"
        )
    )

    container = st.container()

    if app.get_mdata() is None:
        container.page_link(
            "pages/tables.py",
            label="Upload data to get started"
        )
        return

    # If the user has selected a view to edit, show the edit menu
    if st.query_params.get("edit-view") is not None:

        run_edit_view(container)

    else:

        # If the views are editable
        if app.get_edit_views_flag():

            # Show the views with edit buttons
            view_editable(container)

        else:
            # Show the views with no edit buttons
            view_non_editable(container)


def run_edit_view(container: DeltaGenerator):

    # The view to edit
    edit_ix = int(st.query_params.get("edit-view"))

    # Get the list of all views defined in the dataset
    views = app.get_views()

    if len(views) < (edit_ix + 1):
        container.error("No views to edit.")

    # Instantiate the view to edit
    view = make_view(
        ix=edit_ix,
        editable=True,
        **views[edit_ix]
    )

    view.get_data()


def view_editable(container: DeltaGenerator):
    """
    Show the views with edit buttons.
    """

    if not app.get_edit_views_flag():
        st.rerun()

    # All of the views defined in the dataset
    mdata_views = make_views(editable=False)

    for ix, view in enumerate(mdata_views):

        # If the settings are editable
        if app.get_edit_views_flag():
            # Show the name of the view
            container.write(f"#### {ix + 1}. {view.name}")

            # Set up a set of buttons to edit the order of the view
            edit_view(view, container, ix, len(mdata_views))

        # Attach the view to the display
        view.attach(container)

    # Let the user add a new view
    button_add_view()


def view_non_editable(container: DeltaGenerator):

    # All of the views defined in the dataset
    mdata_views = make_views(editable=False)

    for view in mdata_views:

        # Attach the view to the display
        view.attach(container)
