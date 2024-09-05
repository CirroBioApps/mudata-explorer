from mudata import MuData
import os
import streamlit as st
from typing import Optional
from mudata_explorer.app.mdata import set_mdata
from mudata_explorer.app.hash import get_dat_hash, set_mdata_hash
from mudata_explorer.helpers.io import hydrate_uns


def readme(
    folder: str,
    filename="README.md"
) -> str:
    with open(
        os.path.join(
            os.path.abspath(
                os.path.dirname(__file__)
            ),
            folder,
            filename
        ),
        "r"
    ) as fh:
        return fh.read()


def load_mdata(mdata: Optional[MuData] = None, id="main"):
    # If no data was read, stop
    if mdata is None:
        return

    # Get the hash of the data
    _, hash, _ = get_dat_hash(mdata)

    hydrate_uns(mdata)
    set_mdata(mdata, full=True, id=id)
    set_mdata_hash(hash)
    st.switch_page("pages/views.py")
