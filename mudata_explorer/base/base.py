from typing import Union
import anndata as ad
import muon as mu
import numpy as np
import streamlit as st
import json


class MuDataAppHelpers:

    @property
    def refresh_ix(self):
        return st.session_state.get("refresh-ix", 0)

    def get_mdata(self) -> Union[None, mu.MuData]:
        mdata = st.session_state.get("mdata", None)
        if mdata is not None:
            assert isinstance(mdata, mu.MuData)
            views = mdata.uns.get("mudata-explorer-views", [])
            if isinstance(views, str):
                views = json.loads(views)
                mdata.uns["mudata-explorer-views"] = views
                self.set_mdata(mdata)
        return mdata

    def set_mdata(self, mdata: mu.MuData):
        assert isinstance(mdata, mu.MuData)
        st.session_state["mdata"] = mdata

    def setup_mdata(self):
        mdata = mu.MuData({
            '_blank': ad.AnnData(
                X=np.array([[]]),
                obs=[],
                var=[]
            )
        })
        mdata.uns["mudata-explorer-views"] = []
        self.set_mdata(mdata)

    def get_views(self):
        return self.get_mdata().uns.get("mudata-explorer-views", [])

    def set_views(self, views):
        mdata = self.get_mdata()
        mdata.uns["mudata-explorer-views"] = views
        self.set_mdata(mdata)
