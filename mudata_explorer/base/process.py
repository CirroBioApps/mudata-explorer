from typing import List
from mudata_explorer.base.base import MuDataAppHelpers
from streamlit.delta_generator import DeltaGenerator


class Process(MuDataAppHelpers):

    type: str
    name: str
    desc: str
    categories: List[str]

    def run(self, container: DeltaGenerator):
        pass
