from typing import List
import streamlit as st
from streamlit.delta_generator import DeltaGenerator
from mudata_explorer.base.base import MuDataAppHelpers
from mudata_explorer.base.view import View
from mudata_explorer.base.add_data import AddData
from mudata_explorer.helpers import all_views
from mudata_explorer.helpers import get_view_by_type
from mudata_explorer.helpers import make_view


class App(MuDataAppHelpers):

    mudata_uploader: DeltaGenerator
    views_empty: DeltaGenerator
    process_empty: DeltaGenerator
    add_data_empty: DeltaGenerator
    save_load_empty: DeltaGenerator

    def __init__(self):
        self.setup_page()
        self.show_views()
        # self.show_process()
        self.show_add_data()
        # self.show_save_load()

    def setup_page(self):

        st.set_page_config(
            page_title='MuData Explorer',
            layout='centered',
            initial_sidebar_state='auto'
        )

        # Set up tabs
        tabs_labels = ["Views", "Process", "Add Data", "Save/Load"]
        tabs: List[DeltaGenerator] = st.tabs(tabs_labels)

        # Set up an empty in each tab
        for tab_label, tab in zip(tabs_labels, tabs):
            setattr(
                self,
                (
                    tab_label
                    .lower()
                    .replace(" ", "_")
                    .replace("/", "_")
                    + "_empty"
                ),
                tab.empty()
            )

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
        self.refresh_ix_increment("views")

        # Set up a container to display the views
        plots_cont = self.views_empty.container()

        for ix, view in enumerate(self.views):

            # Show the name of the view
            plots_cont.write(f"#### {ix + 1}. {view.name}")

            # Show the display
            view.display(plots_cont)

            # Set up an expander element
            expander = plots_cont.expander("Edit Settings")

            # Populate the form with the view params
            view.inputs(expander)

            # Make a set of columns
            cols = expander.columns([1, 1, 1])

            # Add buttons for editing, reordering, and deleting
            if ix > 0:
                cols[0].button(
                    "Move Up",
                    key=f"up-view-{ix}-{self.refresh_ix('views')}",
                    use_container_width=True,
                    on_click=self.move_up,
                    args=(ix,)
                )
            cols[1].button(
                "Delete",
                key=f"delete-view-{ix}-{self.refresh_ix('views')}",
                use_container_width=True,
                on_click=self.delete_view,
                args=(ix,)
            )
            if ix < (len(self.views) - 1):
                cols[2].button(
                    "Move Down",
                    key=f"down-view-{ix}-{self.refresh_ix('views')}",
                    use_container_width=True,
                    on_click=self.move_down,
                    args=(ix,)
                )

            # Show a horizontal rule
            plots_cont.markdown("---")

        # Let the user add a new view
        self.button_add_view(plots_cont)

    def show_add_data(self):

        # Set up a container to display the views
        add_data_cont = self.add_data_empty.container()

        # Set up the add data object
        add_data = AddData()

        # Show the form
        add_data.show(add_data_cont.empty())

        # Process the inputs
        add_data.process(add_data_cont.empty())

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
                key=f"add-view-button-{ix}-{self.refresh_ix('views')}",
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
