import anndata as ad
import hashlib
import json
import muon as mu
import numpy as np
import streamlit as st
from streamlit.delta_generator import DeltaGenerator
from typing import List, Union
from mudata_explorer.helpers import get_view_by_type


def setup_pages():
    st.set_page_config("MuData Explorer", layout="centered")
    st.sidebar.title("MuData Explorer")
    st.sidebar.page_link("pages/views.py", label="Views")
    st.sidebar.page_link("pages/processes.py", label="Process")
    st.sidebar.page_link("pages/add_data.py", label="Add Data")
    st.sidebar.page_link("pages/save_load.py", label="Save/Load")
    st.sidebar.page_link("pages/settings.py", label="Settings")


def hash_dat(dat):
    """Compute the hash of an object."""
    hash = hashlib.sha256()
    hash.update(dat)
    return hash.hexdigest()


def get_mdata() -> Union[None, mu.MuData]:
    mdata = st.session_state.get("mdata", None)
    if mdata is not None:
        assert isinstance(mdata, mu.MuData)
        views = mdata.uns.get("mudata-explorer-views", [])
        if isinstance(views, str):
            views = json.loads(views)
            mdata.uns["mudata-explorer-views"] = views
            set_mdata(mdata)
    return mdata


def set_mdata(mdata: mu.MuData):
    assert isinstance(mdata, mu.MuData)
    st.session_state["mdata"] = mdata


def setup_mdata():
    mdata = mu.MuData({
        '_blank': ad.AnnData(
            X=np.array([[]]),
            obs=[],
            var=[]
        )
    })
    mdata.uns["mudata-explorer-views"] = []
    set_mdata(mdata)


def summarize_mdata(container: DeltaGenerator):
    mdata = get_mdata()
    if mdata is None:
        return
    if len(mdata.mod) == 1 and "_blank" in mdata.mod:
        return
    container.write("**Current MuData**")
    for mod_name, mod in mdata.mod.items():
        shape = mod.to_df().shape
        container.write(f" - {mod_name} ({shape[0]:,} observations x {shape[1]:,} features)") # noqa


def delete_view(ix: int):
    views = get_views()
    views.pop(ix)
    set_views(views)


def add_view(view_type: str):
    if get_mdata() is None:
        setup_mdata()
    views = get_views()
    views.append(
        get_view_by_type(view_type).template()
    )
    set_views(views)


def get_views() -> List[dict]:
    return get_mdata().uns.get("mudata-explorer-views", [])


def set_views(views):
    mdata = get_mdata()
    mdata.uns["mudata-explorer-views"] = views
    set_mdata(mdata)


def get_settings() -> dict:
    return get_mdata().uns.get("mudata-explorer-settings", {})


def set_settings(settings):
    mdata = get_mdata()
    mdata.uns["mudata-explorer-settings"] = settings
    set_mdata(mdata)
