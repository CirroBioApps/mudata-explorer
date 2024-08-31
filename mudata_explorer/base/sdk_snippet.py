import json
from mudata_explorer.helpers.views import get_views
from mudata_explorer.helpers.params import nest_params
import streamlit as st


@st.dialog("Figure Parameters", width='large')
def show_view_sdk_snippet(ix: int):
    view = get_views()[ix]
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
    return f"""view.{view['type'].replace('-', '_')}(
    mdata,
    **{params_str}
)
"""
