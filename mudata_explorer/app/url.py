import streamlit as st
import requests
from tempfile import NamedTemporaryFile
import muon as mu
from mudata_explorer.app.mdata import set_mdata
from mudata_explorer.helpers.io import hydrate_uns


@st.cache_resource
def _load_url(url: str):
    res = requests.get(url)
    with NamedTemporaryFile(suffix=".h5mu", delete=True) as tmp:
        with open(tmp.file.name, "wb") as f:
            f.write(res.content)
        return mu.read_h5mu(tmp.file.name)


def load_url(url: str, id="main"):
    mdata = _load_url(url)
    hydrate_uns(mdata)
    set_mdata(mdata.copy(), id=id, full=True)
