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

    def __init__(self, ix: int, type: str, **params):
        self.ix = ix
        self.type = type
        self.params = params

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
