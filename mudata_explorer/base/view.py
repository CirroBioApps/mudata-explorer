import streamlit as st
from typing import List
from mudata_explorer.base.base import MuDataAppHelpers
from streamlit.delta_generator import DeltaGenerator


class View(MuDataAppHelpers):

    ix: int
    type: str
    params: dict
    name: str
    desc: str
    categories: List[str]
    logs: List[str] = []
    params: dict = {}
    defaults: dict = {}

    def __init__(
        self,
        ix: int,
        type: str,
        name: str,
        desc: str,
        logs: List[str],
        params: dict,
        on_change: callable,
        refresh: callable
    ):
        self.ix = ix
        self.type = type
        self.name = name
        self.desc = desc
        self.logs = logs
        self.params = {
            kw: params.get(kw, val)
            for kw, val in self.defaults.items()
        }
        self.on_change = on_change
        self.refresh = refresh

    @classmethod
    def template(cls):
        return dict(
            type=cls.type,
            name=cls.name,
            desc=cls.desc,
            logs=cls.logs,
            params=cls.params
        )

    def display(self, container: DeltaGenerator):
        """Primary method which is executed to render the view."""
        pass

    def inputs(self, form: DeltaGenerator):
        """Render any input elements needed for this view within a form."""
        pass

    def process(self):
        pass

    @property
    def refresh_ix(self):
        super().refresh_ix("view")

    @property
    def key_prefix(self):
        return f"view-{self.ix}-{self.refresh_ix}-"

    def param_key(self, key: str):
        return f"{self.key_prefix}{key}"

    def param_kwargs(self, key: str, incl_value=True):
        kwargs = dict(
            key=self.param_key(key),
            on_change=self.on_change,
            args=(self, key,)
        )
        if incl_value:
            kwargs["value"] = self.params.get(key, self.defaults[key])
        return kwargs

    def clear_params(self):
        """Delete any param from the session state which would conflict."""
        to_delete = [
            key
            for key in st.session_state.keys()
            if key.startswith(self.key_prefix)
        ]
        for key in to_delete:
            del st.session_state[key]

    def get_params(self):

        return {
            kw: st.session_state.get(self.param_key(kw), val)
            for kw, val in self.defaults.items()
        }
