from typing import List
import streamlit as st
import muon as mu
from streamlit.delta_generator import DeltaGenerator


class View:

    ix: int
    type: str
    params: dict
    name: str
    desc: str
    categories: List[str]
    processed: bool = False
    logs: List[str] = []
    params: dict = {}

    def __init__(
        self,
        ix: int,
        type: str,
        name: str,
        desc: str,
        processed: bool,
        logs: List[str],
        params: dict
    ):
        self.ix = ix
        self.type = type
        self.name = name
        self.desc = desc
        self.processed = processed
        self.logs = logs
        self.params = params

    @classmethod
    def template(cls):
        return dict(
            type=cls.type,
            name=cls.name,
            desc=cls.desc,
            processed=cls.processed,
            logs=cls.logs,
            params=cls.params
        )

    def inputs(self, container: DeltaGenerator):
        """Render any input elements needed for this view."""
        pass

    def process(self):
        pass

    def display(self, container: DeltaGenerator):
        pass

    @property
    def mdata(self):
        return st.session_state.get("mdata", None)

    @mdata.setter
    def mdata(self, mdata: mu.MuData):
        assert isinstance(mdata, mu.MuData)
        st.session_state["mdata"] = mdata
