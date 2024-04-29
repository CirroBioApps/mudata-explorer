import streamlit as st
from muon import MuData
from typing import List
from mudata_explorer.base.base import MuDataAppHelpers
from streamlit.delta_generator import DeltaGenerator


class Process(MuDataAppHelpers):

    type: str
    name: str
    desc: str
    categories: List[str]
    logs: List[str] = []

    def __init__(
        self,
        on_change: callable
    ):
        self.on_change = on_change

    @classmethod
    def template(cls):
        return dict(
            type=cls.type,
            name=cls.name,
            desc=cls.desc,
        )

    def run(self, container: DeltaGenerator, mdata: MuData):
        pass

    @property
    def refresh_ix(self):
        super().refresh_ix("process")

    @property
    def key_prefix(self):
        return f"process-{self.refresh_ix}-"

    def param_key(self, key: str):
        return f"{self.key_prefix}{key}"
