from mudata_explorer.app.mdata import get_mdata_exists, get_view
from mudata_explorer.app.query_params import get_edit_view_flag
from mudata_explorer.app.session_state import set_show_sidebar_flag
from mudata_explorer.app.sidebar import setup_sidebar
from mudata_explorer.base.sdk_snippet import show_view_sdk_snippet
from mudata_explorer.helpers.assets import all_views, make_view
from mudata_explorer.helpers.assets import asset_categories, asset_dataframe
from mudata_explorer.helpers.assets import filter_by_category
import streamlit as st


def run_edit_view():

    # The view to edit
    edit_ix = get_edit_view_flag()

    # Get the view to edit
    view = get_view(edit_ix)

    # Instantiate the view to edit
    view = make_view(
        ix=edit_ix,
        type=view['type'],
        params=view['params']
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

    # If the selected type is not in the selected category
    if view.type not in df["type"].tolist():

        # Re-instantiate the view to edit
        view = make_view(
            ix=edit_ix,
            type=df["type"].values[0],
            params=view.form.dehydrate()
        )

        # Update the type of the view in the state
        view.save_changes()
        
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
            params=view.form.dehydrate()
        )
        # Update the type of the view in the state
        view.save_changes()
        
    # Render the input form and save any changes in the state
    view.display_form()

    # Show the display
    try:
        view.display()
    except Exception as e:
        st.exception(e)

    st.button(
        ":information_source: SDK Snippet",
        on_click=show_view_sdk_snippet,
        key=f"show-view-sdk-snippet-{view.ix}",
        args=(view.ix,),
        help="Show an example for configuration via SDK."
    )

    if st.button(
        ":page_facing_up: Save Changes",
        key=f"save-changes-{view.ix}"
    ):
        st.switch_page("pages/view_all.py")


def run():

    if not get_mdata_exists():
        st.switch_page("pages/load.py")

    set_show_sidebar_flag(False)

    setup_sidebar("view_details")

    run_edit_view()


if __name__ == "__main__":
    run()
