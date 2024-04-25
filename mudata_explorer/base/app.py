import streamlit as st
from streamlit.delta_generator import DeltaGenerator
from mudata_explorer.base.base import MuDataAppHelpers
from mudata_explorer.base.view import View
from mudata_explorer.helpers import all_views
from mudata_explorer.helpers import get_view_by_type
from mudata_explorer.helpers import make_view


class App(MuDataAppHelpers):

    mudata_uploader: DeltaGenerator

    def __init__(self):
        self.setup_page()
        self.show_views()

    def setup_page(self):

        st.set_page_config(
            page_title='MuData Explorer',
            layout='centered',
            initial_sidebar_state='auto'
        )

        # Set up an empty container to display the views
        self.view_empty = st.empty()

    @property
    def views(self):
        if self.get_mdata() is None:
            return []
        views = self.get_mdata().uns.get("mudata-explorer-views", [])
        assert isinstance(views, list)
        return [
            make_view(
                ix=ix,
                on_change=self.update_view,
                refresh=self.show_views,
                **view
            )
            for ix, view in enumerate(views)
        ]

    def show_views(self):

        # Increment the refresh index
        st.session_state["refresh-ix"] = self.refresh_ix + 1

        # Set up a container to display the views
        view_cont = self.view_empty.container()

        for ix, view in enumerate(self.views):

            # Show the name of the view
            view_cont.write(f"#### {ix + 1}. {view.name}")

            # Show the display
            view.display(view_cont)

            # Set up an expander element
            expander = view_cont.expander("Edit Settings")

            # Populate the form with the view params
            view.inputs(expander)

            # Make a set of columns
            cols = expander.columns([1, 1, 1])

            # Add buttons for editing, reordering, and deleting
            if ix > 0:
                cols[0].button(
                    "Move Up",
                    key=f"up-view-{ix}-{self.refresh_ix}",
                    use_container_width=True,
                    on_click=self.move_up,
                    args=(ix,)
                )
            cols[1].button(
                "Delete",
                key=f"delete-view-{ix}-{self.refresh_ix}",
                use_container_width=True,
                on_click=self.delete_view,
                args=(ix,)
            )
            if ix < (len(self.views) - 1):
                cols[2].button(
                    "Move Down",
                    key=f"down-view-{ix}-{self.refresh_ix}",
                    use_container_width=True,
                    on_click=self.move_down,
                    args=(ix,)
                )

            # Show a horizontal rule
            view_cont.markdown("---")

        # Let the user add a new view
        self.button_add_view(view_cont)

    def update_view(self, view: View, key: str):
        # Get the new value
        val = st.session_state[view.param_key(key)]

        # If it does not match the old value
        if val != view.params[key]:
            # Update the view
            view.params[key] = val
            mdata = self.get_mdata()
            mdata.uns["mudata-explorer-views"][view.ix]["params"] = view.params
            self.set_mdata(mdata)
            self.show_views()

    def move_up(self, ix: int):
        views = self.get_views()
        (
            views[ix - 1],
            views[ix]
        ) = (
            views[ix],
            views[ix - 1]
        )
        self.set_views(views)
        self.show_views()

    def move_down(self, ix: int):
        views = self.get_views()
        (
            views[ix + 1],
            views[ix]
        ) = (
            views[ix],
            views[ix + 1]
        )
        self.set_views(views)
        self.show_views()

    def delete_view(self, ix: int):
        views = self.get_views()
        views.pop(ix)
        self.set_views(views)
        self.show_views()

    def button_add_view(self, container: DeltaGenerator):
        container.write("#### Available Options")
        for ix, view in enumerate(all_views):
            cols = container.columns([1, 5])
            cols[1].write(f"**{view.name}**\n\n{view.desc}")
            cols[0].button(
                "Add",
                help=view.desc,
                key=f"add-view-button-{ix}-{self.refresh_ix}",
                on_click=self.add_view,
                args=(view.type,),
                use_container_width=True
            )

    def add_view(self, view_type: str):
        if self.get_mdata() is None:
            self.setup_mdata()
        views = self.get_views()
        views.append(
            get_view_by_type(view_type).template()
        )
        self.set_views(views)
        self.show_views()
