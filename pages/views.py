from mudata_explorer import app
from mudata_explorer.helpers import all_views, make_view
from mudata_explorer.helpers import asset_categories, asset_type_desc_lists
from mudata_explorer.helpers import filter_by_category
import streamlit as st
from streamlit.delta_generator import DeltaGenerator


def make_views(editable: bool):
    if app.get_mdata() is None:
        return []
    views = app.get_mdata().uns.get("mudata-explorer-views", [])
    assert isinstance(views, list)
    return [
        make_view(
            ix=ix,
            editable=editable,
            **view
        )
        for ix, view in enumerate(views)
    ]


def edit_view(container: DeltaGenerator, ix: int, n_views: int):

    # Make a set of columns
    cols = container.columns([1, 1, 1])

    # Add buttons for editing, reordering, and deleting
    if ix > 0:
        cols[0].button(
            "Move Up",
            use_container_width=True,
            on_click=move_up,
            key=f"move-up-{ix}",
            args=(ix,)
        )
    cols[1].button(
        "Delete",
        use_container_width=True,
        on_click=app.delete_view,
        key=f"delete-view-{ix}",
        args=(ix,)
    )
    if ix < (n_views - 1):
        cols[2].button(
            "Move Down",
            use_container_width=True,
            on_click=move_down,
            key=f"move-down-{ix}",
            args=(ix,)
        )


def button_add_view(container: DeltaGenerator):

    container.write("#### Add a new view")

    # Let the user select the type of view to add
    all_categories = asset_categories(all_views)

    selected_category = container.selectbox(
        "Select a category",
        all_categories
    )

    # Let the user select which view to add, filtering by category
    filtered_views = filter_by_category(all_views, selected_category)

    # Get the assets needed to select from the filtered views
    type_list, desc_list = asset_type_desc_lists(filtered_views)

    selected_desc = container.selectbox(
        "Select a view to add",
        desc_list
    )
    selected_type = type_list[desc_list.index(selected_desc)]

    # Instantiate the selected view type if the user clicks a button
    container.button(
        f"Add {selected_desc}",
        on_click=app.add_view,
        args=(selected_type,),
        use_container_width=True
    )


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


if __name__ == "__main__":

    app.setup_pages()

    container = st.container()

    settings = app.get_settings()

    mdata_views = make_views(editable=settings["editable"])

    for ix, view in enumerate(mdata_views):

        # Show the name of the view
        if settings["editable"]:
            container.write(f"#### {ix + 1}. {view.name}")

        # Attach the view to the display
        view.attach(container)

        # Set up a set of buttons to edit the order of the view
        if settings["editable"]:
            edit_view(view.inputs_container, ix, len(mdata_views))

            # Show a horizontal rule
            container.markdown("---")

    if settings["editable"]:
        # Let the user add a new view
        button_add_view(container)
