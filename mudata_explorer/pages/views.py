from mudata_explorer.helpers.views import get_views, set_views, delete_view, duplicate_view
from mudata_explorer.helpers.add_view import add_view
from mudata_explorer.app.sidebar import get_editable_flag, setup_sidebar
from mudata_explorer.app.mdata import get_mdata
from mudata_explorer.app.query_params import set_edit_view_flag
from mudata_explorer.app.query_params import get_edit_view_flag
from mudata_explorer.base.view import View
from mudata_explorer.base.sdk_snippet import show_view_sdk_snippet
from mudata_explorer.helpers.assets import all_views, make_view
from mudata_explorer.helpers.assets import asset_categories, asset_dataframe
from mudata_explorer.helpers.assets import filter_by_category
import streamlit as st
from streamlit.delta_generator import DeltaGenerator


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


def edit_view(view: View, container: DeltaGenerator, ix: int, n_views: int):

    # If the views are not editable, don't show the edit buttons
    if not get_editable_flag():
        # Instead, just set up the params for the view
        view.display_form()
        return

    # Set up an expand element for all of the rearranging options
    expander = container.expander("More Options", expanded=True)

    # The first button allows the user to provide inputs
    expander.button(
        ":pencil: Edit All",
        help="Edit the inputs for this view.",
        key=f"edit-inputs-{ix}",
        on_click=set_edit_view_flag,
        args=(ix,)
    )

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

    expander.button(
        ":heavy_plus_sign: Duplicate",
        on_click=duplicate_view,
        key=f"duplicate-view-{ix}",
        args=(ix,),
        help="Make a copy of this view."
    )

    if ix < (n_views - 1):
        expander.button(
            ":arrow_down_small: Move Down",
            on_click=move_down,
            key=f"move-down-{ix}",
            args=(ix,),
            help="Move this view down in the list."
        )

    # Let the user add text or a figure above
    expander.button(
        f":pencil2: Add Text Above",
        on_click=button_add_view_callback,
        key=f"add-text-above-{ix}",
        args=("markdown",),
        kwargs=dict(ix=ix)
    )

    expander.button(
        f":chart_with_upwards_trend: Add Figure Above",
        on_click=button_add_view_callback,
        key=f"add-figure-above-{ix}",
        args=("plotly-scatter",),
        kwargs=dict(ix=ix)
    )


def button_add_view():

    # Let the user add a text block
    st.button(
        f"Add Text",
        on_click=button_add_view_callback,
        args=("markdown",),
        kwargs=dict(ix=-1)
    )
    # Let the user add a scatter plot
    st.button(
        f"Add Figure",
        on_click=button_add_view_callback,
        args=("plotly-scatter",),
        kwargs=dict(ix=-1)
    )


def button_add_view_callback(selected_type: str, ix=-1):
    add_view(selected_type, ix=ix)
    if ix == -1:
        set_edit_view_flag(len(get_views()) - 1)
    else:
        set_edit_view_flag(ix)


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

    setup_sidebar(
        edit_views=True,
        page_layout=(
            "centered"
            if get_edit_view_flag() is None
            and get_editable_flag() is False
            else "wide"
        )
    )

    if get_mdata() is None:
        st.page_link(
            "pages/tables.py",
            label="Upload data to get started",
            icon=":material/table:"
        )
        return

    # If the user has selected a view to edit, show the edit menu
    if get_edit_view_flag() is not None:

        run_edit_view()

    else:

        # If the views are editable
        if get_editable_flag():

            # Show the views with edit buttons
            view_editable()

        else:
            # Show the views with no edit buttons
            view_non_editable()


def run_edit_view():

    # The view to edit
    edit_ix = get_edit_view_flag()

    # Get the list of all views defined in the dataset
    views = get_views()

    if len(views) < (edit_ix + 1):
        st.error("No views to edit.")

    # Instantiate the view to edit
    view = make_view(
        ix=edit_ix,
        type=views[edit_ix]['type'],
        params=views[edit_ix]['params']
    )

    # Let the user modify the view type
    all_categories = asset_categories(all_views)

    st.markdown("#### Figure Type")
    selected_category = st.selectbox(
        "Category",
        all_categories,
        index=all_categories.index(view.category)
    )

    # Let the user select which view to add, filtering by category
    filtered_views = filter_by_category(all_views, selected_category)

    # Get the assets needed to select from the filtered views
    df = asset_dataframe(filtered_views)

    selected_name = st.selectbox(
        "Figure Type",
        df["name"].tolist(),
        index=df["name"].tolist().index(view.name)
    )
    selected_ix = df["name"].tolist().index(selected_name)
    selected_type = df["type"].tolist()[selected_ix]
    help_text = df["help_text"].tolist()[selected_ix]

    if help_text is not None:
        st.markdown(help_text)

    # If the selected type is different
    if selected_type != view.type:
        
        # Re-instantiate the view to edit
        view = make_view(
            ix=edit_ix,
            type=selected_type,
            params=views[edit_ix]['params']
        )

    # Render the input form
    view.display_form()

    # Show the display
    view.display()

    st.button(
        ":information_source: SDK Snippet",
        on_click=show_view_sdk_snippet,
        key=f"show-view-sdk-snippet-{view.ix}",
        args=(view.ix,),
        help="Show an example for configuration via SDK."
    )

    st.button(
        ":page_facing_up: Save Changes",
        key=f"save-changes-{view.ix}",
        on_click=view.save_changes_button
    )


def view_editable():
    """
    Show the views with edit buttons.
    """

    if not get_editable_flag():
        st.rerun()

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
            view.display_form()

        # Let the user run the method, catching any errors
        if not view.form.complete:
            display.write("Please complete all input fields")
            continue

        # Set up a set of buttons to edit the order of the view
        edit_view(view, controls, ix, len(mdata_views))

        # Now make the display, catching any errors
        try:
            with display:
                view.display()
        except Exception as e:
            # Log the full traceback of the exception
            display.exception(e)

    if len(mdata_views) > 0:
        st.write("---")

    # Let the user add a new view
    button_add_view()


def view_non_editable():

    # All of the views defined in the dataset
    mdata_views = make_views()

    for view in mdata_views:

        # Attach the view to the display
        view.attach()
