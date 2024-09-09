import json
from mudata_explorer.helpers.params import nest_params
from mudata_explorer.app import mdata as mdata_funcs
from mudata_explorer.helpers.io import json_safe
import mudata as mu
import streamlit as st


def update_process_on_change(kw) -> None:
    val = st.session_state[f"process-{kw}"]
    update_process(kw, val)


def update_process(kw, val, id="main") -> None:
    process = mdata_funcs.get_process(id=id)
    process[kw] = val
    mdata_funcs.set_process(process, id=id)


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


