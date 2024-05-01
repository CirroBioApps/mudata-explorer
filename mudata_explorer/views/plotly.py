from typing import Tuple, Union
import anndata as ad
from mudata_explorer import app
from mudata_explorer.base.view import View
import pandas as pd
import plotly.express as px
from streamlit.delta_generator import DeltaGenerator


class Plotly(View):

    categories = ["Plotting"]
    defaults = {
        "modality": None,
        "x": None,
        "y": None,
        "z": None,
        "sort_by": None,
        "query": "",
        "use_color": False,
        "color": None,
        "color_continuous_scale": "blues"
    }

    def get_data(
        self,
        container: DeltaGenerator,
        keys=["x", "y", "color"]
    ) -> Union[None, pd.DataFrame]:

        for kw in keys:
            assert kw in self.params, f"Missing parameter: {kw}"

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

        # Make a single DataFrame with all of the data for this modality
        df = app.make_modality_df(mdata, self.params["modality"])

        # Let the user optionally filter samples
        if settings["editable"]:
            self.params["query"] = container.text_input(
                "Filter samples (optional)",
                help="Enter a query to filter samples (using metadata or data).",
                **self.input_value_kwargs("query")
            )
        if self.params["query"] is not None and len(self.params["query"]) > 0:
            df = df.query(self.params["query"])

            if settings["editable"]:
                container.write(
                    f"Filtered data: {df.shape[0]:,} rows x {df.shape[1]:,} columns."
                )

        # If there is no data, return
        if df.shape[0] == 0 or df.shape[1] == 0:
            container.write(f"No data available for {self.params['modality']}.")
            return

        cols = list(df.columns.values)

        all_params = [
            dict(kw="x", label="Select x-axis"),
            dict(kw="y", label="Select y-axis"),
            dict(kw="z", label="Select z-axis"),
            dict(kw="sort_by", label="Sort data by")
        ]

        choice_ix = 0
        for param in all_params:
            kw = param["kw"]
            if self.params.get(kw) is None and kw in keys:
                self.update_view_param(kw, cols[choice_ix])
                if choice_ix + 1 < len(cols):
                    choice_ix += 1

            if settings["editable"] and kw in keys:
                self.params[kw] = container.selectbox(
                    param["label"],
                    cols,
                    **self.input_selectbox_kwargs(kw, cols)
                )

        if settings["editable"] and "color" in keys:
            if container.checkbox(
                "Use color",
                **self.input_value_kwargs("use_color")
            ):
                self.params["color"] = container.selectbox(
                    "Select color",
                    cols,
                    **self.input_selectbox_kwargs("color", cols)
                )
                self.params["color_continuous_scale"] = container.selectbox(
                    "Select color scale",
                    px.colors.named_colorscales(),
                    **self.input_selectbox_kwargs(
                        "color_continuous_scale",
                        px.colors.named_colorscales()
                    )
                )
            else:
                self.update_view_param("color", None)

        return df


class PlotlyScatter(Plotly):

    type = "plotly-scatter"
    name = "Scatterplot (Plotly)"
    desc = "Display a two dimensional distribution of data using Plotly."

    def display(self, container: DeltaGenerator):

        data = self.get_data(
            container,
            keys=["x", "y", "color"]
        )
        if data is None:
            return

        cols = [self.params["x"], self.params["y"]]
        if self.params["color"] is not None:
            cols.append(self.params["color"])

        cols = list(set(cols))

        data = data.reindex(columns=cols).dropna()

        fig = px.scatter(
            data,
            x=self.params["x"],
            y=self.params["y"],
            color=self.params["color"],
            color_continuous_scale=self.params["color_continuous_scale"]
        )

        container.plotly_chart(fig)


class PlotlyScatter3D(Plotly):

    type = "plotly-scatter-3d"
    name = "Scatterplot 3D (Plotly)"
    desc = "Display a three dimensional distribution of data using Plotly."

    def display(self, container: DeltaGenerator):

        data = self.get_data(
            container,
            keys=["x", "y", "z", "color"]
        )
        if data is None:
            return

        cols = [self.params["x"], self.params["y"], self.params["z"]]
        if self.params["color"] is not None:
            cols.append(self.params["color"])

        cols = list(set(cols))

        data = data.reindex(columns=cols).dropna()

        fig = px.scatter_3d(
            data,
            x=self.params["x"],
            y=self.params["y"],
            z=self.params["z"],
            color=self.params["color"],
            color_continuous_scale=self.params["color_continuous_scale"]
        )

        container.plotly_chart(fig)


class PlotlyLine(Plotly):

    type = "plotly-line"
    name = "Line Plot (Plotly)"
    desc = "Display a series of data as a line graph using Plotly."

    def display(self, container: DeltaGenerator):

        data = self.get_data(
            container,
            keys=["x", "y", "sort_by", "color"]
        )
        if data is None:
            return

        cols = [self.params["x"], self.params["y"], self.params["sort_by"]]
        if self.params["color"] is not None:
            cols.append(self.params["color"])

        cols = list(set(cols))

        data = data.reindex(columns=cols).dropna()
        data = data.sort_values(self.params["sort_by"])

        fig = px.line(
            data,
            x=self.params["x"],
            y=self.params["y"],
            color=self.params["color"]
        )

        container.plotly_chart(fig)
