from mudata_explorer.app.mdata import get_view, get_views, set_views, get_mdata_exists
from mudata_explorer.app.query_params import set_edit_view_flag
from mudata_explorer.app.session_state import set_show_sidebar_flag
from mudata_explorer.app.write_image import dialog_write_image
from mudata_explorer.app.sidebar import setup_sidebar
from mudata_explorer.helpers.add_view import add_view
from mudata_explorer.helpers.assets import make_view
from mudata_explorer.helpers.views import delete_view, duplicate_view
from streamlit.delta_generator import DeltaGenerator
import streamlit as st


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


def view_editable():
    """
    Show the views with edit buttons.
    """

    # All of the views defined in the dataset
    mdata_views = make_views()

    for ix, view in enumerate(mdata_views):
        if ix > 0:
            st.write("---")

        # Set up two columns
        # The controls will be shown on the left,
        # and the display on the right
        controls, display = st.columns([1, 2])

        # Show the name of the view
        controls.write(f"#### {ix + 1}. {view.name}")

        # In the controls columns
        with controls:
            # Expose any params which can be configured in the sidebar
            view.runtime_options()
            # Show any sidebar parameters in the controls container
            try:
                view.display_form()
            except Exception as e:
                st.exception(e)

        # Let the user run the method, catching any errors
        if not view.form.complete:

            display.write("Please complete all input fields")

        else:

            # Now make the display, catching any errors
            try:
                with display:
                    view.display()
            except Exception as e:
                # Log the full traceback of the exception
                display.exception(e)

        # Set up a set of buttons to edit the order of the view
        edit_view(controls, ix, len(mdata_views))

    if len(mdata_views) > 0:
        st.write("---")

    # Let the user add a new view
    button_add_view()


def edit_view(container: DeltaGenerator, ix: int, n_views: int):

    # Set up an expand element for all of the rearranging options
    expander = container.expander("More Options", expanded=True)

    # The first button allows the user to provide inputs
    if expander.button(
        ":pencil: Edit All",
        help="Edit the inputs for this view.",
        key=f"edit-inputs-{ix}"
    ):
        set_edit_view_flag(ix)

    # For every view past the first
    if ix > 0:
        expander.button(
            ":arrow_up_small: Move Up",
            on_click=move_up,
            key=f"move-up-{ix}",
            args=(ix,),
            help="Move this view up in the list.",
            disabled=ix == 0
        )

    if expander.button(
        ":wastebasket: Delete",
        key=f"delete-view-{ix}",
        args=(ix,),
        help="Remove this view."
    ):
        delete_view(ix)
        st.rerun()

    if expander.button(
        ":heavy_plus_sign: Duplicate",
        key=f"duplicate-view-{ix}",
        help="Make a copy of this view."
    ):
        duplicate_view(ix)
        set_edit_view_flag(ix+1)

    if ix < (n_views - 1):
        expander.button(
            ":arrow_down_small: Move Down",
            on_click=move_down,
            key=f"move-down-{ix}",
            args=(ix,),
            help="Move this view down in the list."
        )

    # Let the user add text or a figure above
    if expander.button(
        f":pencil2: Add Text Above",
        key=f"add-text-above-{ix}"
    ):
        button_add_view_callback("markdown", ix=ix)

    if expander.button(
        f":chart_with_upwards_trend: Add Figure Above",
        key=f"add-figure-above-{ix}"
    ):
        button_add_view_callback("plotly-scatter", ix=ix)

    # Optionally show a button used to write the image to a file
    with expander:
        button_write_image(ix)


def button_add_view():

    # Let the user add a text block
    if st.button(f"Add Text"):
        button_add_view_callback("markdown", ix=-1)

    # Let the user add a scatter plot
    if st.button(f"Add Figure"):
        button_add_view_callback("plotly-scatter", ix=-1)


def button_add_view_callback(selected_type: str, ix=-1):
    add_view(selected_type, ix=ix)
    if ix == -1:
        set_edit_view_flag(len(get_views()) - 1)
    else:
        set_edit_view_flag(ix)


def button_write_image(ix: int, id="main"):

    # Get the information on the view
    view = get_view(ix=ix, id=id)

    # Instantiate the View object
    view = make_view(
        ix=ix,
        type=view['type'],
        params=view['params']
    )
    view.params = view.form.dump()

    # If the view has a write_image method
    if not hasattr(view, "write_image"):
        return

    st.button(
        ":camera: Save Image",
        key=f"save-image-{ix}",
        on_click=dialog_write_image,
        args=(view,)
    )


def move_up(ix: int):
    views = get_views()
    (
        views[ix - 1],
        views[ix]
    ) = (
        views[ix],
        views[ix - 1]
    )
    set_views(views)


def move_down(ix: int):
    views = get_views()
    (
        views[ix + 1],
        views[ix]
    ) = (
        views[ix],
        views[ix + 1]
    )
    set_views(views)


def run():

    if not get_mdata_exists():
        st.switch_page("pages/load.py")

    set_show_sidebar_flag(True)

    setup_sidebar("view_sidebar")

    # Show the views with edit buttons
    view_editable()


if __name__ == "__main__":
    run()
