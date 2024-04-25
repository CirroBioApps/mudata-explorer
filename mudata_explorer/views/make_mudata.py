import anndata as ad
from io import StringIO
import muon as mu
import pandas as pd
from mudata_explorer.base.view import View
from streamlit.delta_generator import DeltaGenerator
import streamlit as st


class MakeMuData(View):

    type = "make-mudata"
    name = "Upload Data Tables"
    desc = "Create a MuData object by uploading a set of CSV/TSV tables."
    categories = ["Data Processing"]
    defaults = {
        "obs_file": None,
        "obs": None,
        "mod_name": "Measurement",
        "mod_file": None,
        "mod": None,
        "mod_dict": dict()
    }

    def display(self, container: DeltaGenerator):
        if self.params["obs"] is None:
            obs = self.get_mdata().obs
        else:
            obs = pd.read_csv(StringIO(self.params["obs"]))

        container.write(f"Metadata: {obs.shape[0]:,} rows x {obs.shape[1]:,} columns.") # noqa

        if self.params["mod"] is None:
            container.write("No measurement data uploaded.")
            return

        mod = pd.read_csv(StringIO(self.params["mod"]))
        container.write(f"Measurement: {mod.shape[0]:,} rows x {mod.shape[1]:,} columns.") # noqa
        container.write(f"Measurement Name: {self.params['mod_name']}")

        # Find the overlapping index labels between the two tables
        overlap = set(obs.index).intersection(set(mod.index))

        container.write(f"No. of overlapping observations: {len(overlap):,}")

        if len(overlap) == 0:
            return

        if container.button("Create MuData", key=self.param_key("create")):
            obs = obs.loc[list(overlap)]
            mod = mod.loc[list(overlap)]

            mdata = mu.MuData({
                self.params['mod_name']: ad.AnnData(
                    obs=obs,
                    X=mod
                )
            })
            mdata.update_obs()
            views = self.get_views()
            # Delete the data for this view
            for kw in ["obs", "mod", "mod_name"]:
                del views[self.ix]["params"][kw]
            # Set the views in the new MuData object
            mdata.uns["mudata-explorer-views"] = views
            self.set_mdata(mdata)
            self.refresh()

    def inputs(self, form: DeltaGenerator):
        form.file_uploader(
            "Observation Metadata (.obs)",
            help="Provide a CSV/TSV where the first column is a unique identifier for each observation.", # noqa
            key=self.param_key("obs_file"),
            on_change=self.read_table,
            args=("obs_file", "obs")
        )

        form.text_input(
            "Measurement Name",
            help="Enter the name of the measurement.",
            **self.param_kwargs("mod_name")
        )

        form.file_uploader(
            "Measurement Data (.X)",
            help="Provide a CSV/TSV where the first column is a unique identifier for each observation.", # noqa
            key=self.param_key("mod_file"),
            on_change=self.read_table,
            args=("mod_file", "mod")
        )

    def read_table(self, source_key, dest_key):
        file = st.session_state.get(self.param_key(source_key), None)
        if file is None:
            return

        if file.name.endswith("xlsx"):
            df = pd.read_excel(file)
        elif file.name.endswith("csv"):
            df = pd.read_csv(file)
        elif file.name.endswith("tsv"):
            df = pd.read_csv(
                file,
                sep="\t"
            )
        else:
            st.error("File must be a CSV or TSV.")
            return

        # The first column must only have unique values
        if not df.iloc[:, 0].is_unique:
            st.error("The first column must have unique values.")
            return

        # Add the data to params
        st.session_state[self.param_key(dest_key)] = df.to_csv(index=None)

        # Update the app
        self.on_change(self, dest_key)
