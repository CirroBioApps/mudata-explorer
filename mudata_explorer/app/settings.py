from mudata_explorer.app import mdata
from mudata_explorer.helpers.io import json_safe, validate_json
from mudata_explorer.app.mdata import set_mdata
import muon as mu
import json
import streamlit as st


def get_settings() -> dict:
    if not mdata.has_mdata():
        settings = {}
    else:
        settings = json_safe(
            mdata.get_mdata().uns.get("mudata-explorer-settings", {})
        )

    return settings


def set_settings(settings: dict):
    mdata = mdata.get_mdata()
    assert isinstance(mdata, mu.MuData), type(mdata)

    # Make sure that the data is JSON serializable
    settings = validate_json(settings)

    if json.dumps(settings) != json.dumps(get_settings()):

        mdata.uns["mudata-explorer-settings"] = settings
        set_mdata(mdata)
        st.experimental_rerun()
