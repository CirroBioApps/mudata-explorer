import json
from mudata_explorer.app.mdata import get_view
from mudata_explorer.helpers.params import nest_params
import streamlit as st


@st.dialog("Figure Parameters", width='large')
def show_view_sdk_snippet(ix: int):
    view = get_view(ix)
    st.code(sdk_snippet(view))


def sdk_snippet(view: dict):
    assert "type" in view.keys()
    assert "params" in view.keys()
    params = nest_params(view["params"])
    params_str = (
        json.dumps(params, indent=4)
        .replace('false', 'False')
        .replace('true', 'True')
        .replace('null', 'None')
        .replace("\n", "\n    ")
    )
    return f"""from mudata_explorer.sdk import view
view.{view['type'].replace('-', '_')}(
    mdata,
    **{params_str}
)
"""
