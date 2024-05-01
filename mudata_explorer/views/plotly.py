from typing import Tuple, Union
import anndata as ad
from mudata_explorer import app
from mudata_explorer.base.view import View
import pandas as pd
import plotly.express as px
from streamlit.delta_generator import DeltaGenerator


class PlotlyScatter(View):

    type = "plotly-scatter"
    name = "Scatterplot (Plotly)"
    desc = "Display a two dimensional distribution of data using Plotly."
    categories = ["Plotting"]
    defaults = {
        "modality": None,
        "slot": "X",
        "x": None,
        "y": None,
        "query": "",
        "use_color": False,
        "color": None
    }

    def list_slots(self, adata: ad.AnnData):
        return (
            ["X"] + 
            ["obsm." + kw for kw in list(adata.obsm.keys())] +
            ["obsp." + kw for kw in list(adata.obsp.keys())]
        )

    def get_slot_data(self, adata: ad.AnnData, slot: str):
        if slot == "X":
            return adata.to_df()
        elif slot.startswith("obsm."):
            return pd.DataFrame(adata.obsm[slot.split(".", 1)[1]])
        elif slot.startswith("obsp."):
            return pd.DataFrame(adata.obsp[slot.split(".", 1)[1]])
        else:
            raise ValueError(f"Unknown slot: {slot}")

    def get_data(
        self,
        container: DeltaGenerator,
        keys=["x", "y", "color"]
    ) -> Union[None, Tuple[pd.DataFrame, dict]]:

        mdata = app.get_mdata()
        settings = app.get_settings()

        if mdata is None or mdata.shape[0] == 0:
            container.write("No MuData object available.")
            return

        # Select the modality to use
        if settings["editable"]:
            modality_options = list(mdata.mod.keys())
            container.selectbox(
                "Select modality",
                modality_options,
                **self.input_selectbox_kwargs("modality", modality_options)
            )

        adata: ad.AnnData = mdata.mod[self.params["modality"]]

        if settings["editable"]:
            slots = self.list_slots(adata)
            container.selectbox(
                "Select slot",
                slots,
                **self.input_selectbox_kwargs("slot", slots)
            )

        df = self.get_slot_data(adata, self.params["slot"])

        # If there is no data, return
        if df.shape[0] == 0 or df.shape[1] == 0:
            container.write(f"No data available for {self.params['modality']}.")
            return

        # Let the user optionally filter samples
        if settings["editable"]:
            container.text_input(
                "Filter samples (optional)",
                help="Enter a query to filter samples (using metadata or data).",
                **self.input_value_kwargs("query")
            )
        if self.params["query"] is not None and len(self.params["query"]) > 0:
            df = (
                df
                .merge(mdata.obs, left_index=True, right_index=True)
                .query(self.params["query"])
                .dropna()
            )
            if settings["editable"]:
                container.write(
                    f"Filtered data: {df.shape[0]:,} rows x {df.shape[1]:,} columns."
                )

        if df.shape[0] == 0:
            return

        # Add the obs metadata to the data available for plotting
        df = df.merge(mdata.obs, left_index=True, right_index=True)

        cols = list(df.columns.values)

        if self.params["x"] is None and "x" in keys:
            self.update_view_param("x", cols[0])
        if self.params["y"] is None and "y" in keys:
            self.update_view_param("y", cols[1])
        if self.params["color"] is None and "color" in keys:
            self.update_view_param("color", cols[2])

        if settings["editable"] and "x" in keys:
            container.selectbox(
                "Select x-axis",
                cols,
                **self.input_selectbox_kwargs("x", cols)
            )

        if settings["editable"] and "y" in keys:
            container.selectbox(
                "Select y-axis",
                df.columns,
                **self.input_selectbox_kwargs("y", cols)
            )

        if settings["editable"] and "color" in keys:
            if container.checkbox(
                "Use color",
                **self.input_value_kwargs("use_color")
            ):
                container.selectbox(
                    "Select color",
                    cols,
                    **self.input_selectbox_kwargs("color", cols)
                )

        kwargs = dict(
            x=self.params["x"],
            y=self.params["y"],
            color=self.params["color"] if self.params["use_color"] else None
        )

        return df, kwargs

    def display(self, container: DeltaGenerator):

        data = self.get_data(
            container,
            keys=["x", "y", "color"]
        )
        if data is None:
            return

        df, kwargs = data

        fig = px.scatter(
            df,
            **kwargs
        )

        container.plotly_chart(fig)
