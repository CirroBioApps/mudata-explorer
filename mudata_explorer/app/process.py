import json
from mudata_explorer.app.mdata import has_mdata, get_mdata, setup_mdata, set_mdata
from mudata_explorer.helpers.io import json_safe
import muon as mu
import streamlit as st


def get_process() -> dict:
    if not has_mdata():
        return {}
    mdata = get_mdata()
    assert isinstance(mdata, mu.MuData), type(mdata)
    return json_safe(mdata.uns.get("mudata-explorer-process", {}))


def set_process(process: dict) -> None:
    if not has_mdata():
        setup_mdata()

    mdata = get_mdata()
    assert mdata is not None
    assert isinstance(mdata, mu.MuData), type(mdata)
    mdata.uns["mudata-explorer-process"] = process
    set_mdata(mdata)


def update_process_on_change(kw) -> None:
    val = st.session_state[f"process-{kw}"]
    update_process(kw, val)


def update_process(kw, val) -> None:
    process = get_process()
    process[kw] = val
    set_process(process)


def process_sdk_snippet(prov: dict):
    assert "process" in prov.keys()
    assert "params" in prov.keys()
    params = nest_params(prov["params"])
    params_str = (
        json.dumps(params, indent=4)
        .replace('false', 'False')
        .replace('true', 'True')
        .replace('null', 'None')
        .replace("\n", "\n    ")
    )
    return f"""process.{prov['process'].replace('-', '_')}(
    mdata,
    **{params_str}
)
"""


def nest_params(params: dict):
    output = dict()
    for key, value in params.items():
        if '.' in key:
            keys = key.split('.', 1)
            if keys[0] not in output:
                output[keys[0]] = dict()
            output[keys[0]][keys[1]] = value
        else:
            output[key] = value
    return {
        key: nest_params(value) if isinstance(value, dict) else value
        for key, value in output.items()
    }

