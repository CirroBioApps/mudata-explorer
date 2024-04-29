import anndata as ad
from io import StringIO
import muon as mu
import pandas as pd
from mudata_explorer.base.base import MuDataAppHelpers
from mudata_explorer.helpers import read_table
from mudata_explorer.helpers import sanitize_types
from streamlit.delta_generator import DeltaGenerator
import streamlit as st


class AddData(MuDataAppHelpers):

    def __init__(self):
        pass

    @property
    def key_prefix(self):
        return f"add-data-{self.refresh_ix('add-data')}-"

    def param_key(self, key: str):
        return f"{self.key_prefix}{key}"

    def add_mudata(self, obs: pd.DataFrame, mod: pd.DataFrame, mod_name: str):
        overlap = set(obs.index).intersection(set(mod.index))
        obs = obs.loc[list(overlap)]
        mod = mod.loc[list(overlap)]

        adata = ad.AnnData(obs=obs, X=mod)

        mdata = self.get_mdata()
        if mdata is None:
            mdata = mu.MuData({
                mod_name: adata
            })

        else:

            mdata = mu.MuData({
                **{
                    kw: adata
                    for kw, adata in mdata.mod.items()
                    if kw != '_blank'
                },
                **{
                    mod_name: adata
                }
            })

        mdata.update_obs()
        self.set_mdata(mdata)

    def show(self, empty: DeltaGenerator):
        container = empty.container()
        self.summarize_mdata(container)

        # Get input from the user
        container.file_uploader(
            "Observation Metadata (.obs)",
            help="Provide a CSV/TSV where the first column is a unique identifier for each observation.", # noqa
            key=self.param_key("obs_file")
        )

        container.text_input(
            "Measurement Name",
            value="Measurement",
            help="Enter the name of the measurement.",
            key=self.param_key("mod_name")
        )

        container.file_uploader(
            "Measurement Data (.X)",
            help="Provide a CSV/TSV where the first column is a unique identifier for each observation.", # noqa
            key=self.param_key("mod_file")
        )

    def process(self, empty: DeltaGenerator):
        container = empty.container()
        obs_file = st.session_state.get(self.param_key("obs_file"))
        obs = read_table(obs_file, container)
        mod_name = st.session_state.get(self.param_key("mod_name"))
        mod_file = st.session_state.get(self.param_key("mod_file"))
        mod = read_table(mod_file, container)

        if obs is None:
            container.write("No metadata uploaded")

        else:
            obs = sanitize_types(obs, container)
            container.write(f"Metadata: {obs.shape[0]:,} rows x {obs.shape[1]:,} columns.") # noqa

        if mod is None:
            container.write("No measurement data uploaded.")

        else:

            mod = sanitize_types(mod, container)
            container.write(f"{mod_name}: {mod.shape[0]:,} rows x {mod.shape[1]:,} columns.") # noqa

        if mod is None or obs is None:
            return

        # Find the overlapping index labels between the two tables
        overlap = set(obs.index).intersection(set(mod.index))

        container.write(f"No. of overlapping observations: {len(overlap):,}")

        if len(overlap) == 0:
            return

        container.button(
            "Add to MuData",
            key=self.param_key("create"),
            on_click=self.add_mudata,
            args=(obs, mod, mod_name)
        )
