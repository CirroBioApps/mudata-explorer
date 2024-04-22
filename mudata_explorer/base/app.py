from typing import Union
import muon as mu
import streamlit as st
from streamlit.delta_generator import DeltaGenerator
from mudata_explorer.base.view import View


class App:

    mudata_uploader: DeltaGenerator
    mdata: Union[None, mu.MuData]

    def __init__(self):
        self.setup_page()
        self.read_inputs()
        self.show_outputs()

    @property
    def mdata(self):
        return st.session_state.get("mdata", None)

    @mdata.setter
    def mdata(self, mdata: mu.MuData):
        assert isinstance(mdata, mu.MuData)
        st.session_state["mdata"] = mdata

    @property
    def views(self):
        if self.mdata is None:
            return []
        views = self.mdata.uns["mudata-explorer-views"]
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

        self.mudata_uploader = st.sidebar.file_uploader(
            label='Upload MuData file',
            type=['h5mu'],
            key="mudata_uploader"
        )

    def read_inputs(self):
        if st.session_state["mudata_uploader"]:
            self.mdata = mu.read_h5mu(
                st.session_state["mudata_uploader"]
            )
        else:
            self.mdata = None

    def show_outputs(self):
        for view in self.views:
            view.display()
