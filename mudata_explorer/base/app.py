import anndata as ad
from typing import Union
import muon as mu
import numpy as np
import streamlit as st
from streamlit.delta_generator import DeltaGenerator
from mudata_explorer.base.view import View
from mudata_explorer.helpers import all_views
from mudata_explorer.helpers import get_view_by_type


class App:

    mudata_uploader: DeltaGenerator
    mdata: Union[None, mu.MuData]

    def __init__(self):
        self.setup_page()
        self.show_views()

    @property
    def mdata(self):
        return st.session_state.get("mdata", None)

    @mdata.setter
    def mdata(self, mdata: mu.MuData):
        assert isinstance(mdata, mu.MuData)
        st.session_state["mdata"] = mdata

    def setup_mdata(self):
        self.mdata = mu.MuData({
            'blank': ad.AnnData(
                X=np.array([[]]),
                obs=[],
                var=[]
            )
        })
        self.mdata.uns["mudata-explorer-views"] = []

    @property
    def views(self):
        if self.mdata is None:
            return []
        views = self.mdata.uns.get("mudata-explorer-views", [])
        assert isinstance(views, list)
        return [
            View(
                ix=ix,
                **view
            )
            for ix, view in enumerate(views)
        ]

    def setup_page(self):

        st.set_page_config(
            page_title='MuData Explorer',
            layout='centered',
            initial_sidebar_state='auto'
        )

        # Set up an empty container to display the views
        self.view_empty = st.empty()

    def show_views(self):

        # Increment the refresh index
        st.session_state["refresh-ix"] = self.refresh_ix + 1

        # Set up a container to display the views
        view_cont = self.view_empty.container()

        view_cont.write(self.views)

        for ix, view in enumerate(self.views):
            # Show the logs
            view_cont.write("\n".join(view.logs))
            # Show the inputs
            view.inputs(view_cont)

            # Show the display
            if view.processed:
                view.display(view_cont)

            # Let the user delete the view
            self.button_delete_view(view_cont, ix)

        # Let the user add a new view
        self.button_add_view(view_cont)

    @property
    def refresh_ix(self):
        return st.session_state.get("refresh-ix", 0)

    def button_delete_view(self, container: DeltaGenerator, ix: int):
        container.button(
            "Delete",
            key=f"delete-view-{ix}-{self.refresh_ix}",
            on_click=self.delete_view,
            args=(ix,)
        )

    def delete_view(self, ix: int):
        mdata = self.mdata
        mdata.uns["mudata-explorer-views"].pop(ix)
        self.mdata = mdata
        self.show_views()

    def button_add_view(self, container: DeltaGenerator):
        container.write("### Add")
        for ix, view in enumerate(all_views):
            container.button(
                view.name,
                help=view.desc,
                key=f"add-view-button-{ix}-{self.refresh_ix}",
                on_click=self.add_view,
                args=(view.type,)
            )

    def add_view(self, view_type: str):
        if self.mdata is None:
            self.setup_mdata()
        mdata = self.mdata
        mdata.uns["mudata-explorer-views"].append(
            get_view_by_type(view_type).template()
        )
        self.mdata = mdata
        self.show_views()
