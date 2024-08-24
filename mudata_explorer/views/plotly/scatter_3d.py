import plotly.express as px
import streamlit as st
from mudata_explorer.views.plotly.base import Plotly


class PlotlyScatter3D(Plotly):

    type = "plotly-scatter-3d"
    name = "Scatterplot 3D (Plotly)"
    help_text = "Display a three dimensional distribution of data using Plotly"
    schema = {
        "data": {
            "type": "dataframe",
            "columns": {
                "x": {"label": "x-axis"},
                "y": {"label": "y-axis"},
                "z": {"label": "z-axis"},
                "size": {"label": "size", "optional": True},
                "color": {
                    "label": "Color",
                    "optional": True,
                    "colorscale": True
                }
            },
            "sidebar": True,
            "query": True,
        },
        "scale_options": {
            "type": "object",
            "label": "Scale Options",
            "properties": {
                "log_x": {
                    "type": "boolean",
                    "label": "Log Scale - X Axis",
                    "sidebar": True
                },
                "log_y": {
                    "type": "boolean",
                    "label": "Log Scale - Y Axis",
                    "sidebar": True
                },
                "log_z": {
                    "type": "boolean",
                    "label": "Log Scale - Z Axis",
                    "sidebar": True
                }
            }
        }
    }

    def display(self):

        data, colorscale = self.fetch_dataframe("data")
        if data is None:
            st.write("Please select a data table")
            return

        fig = px.scatter_3d(
            data.reset_index(),
            x="x",
            y="y",
            z="z",
            hover_name=(
                "index" if data.index.name is None else
                data.index.name
            ),
            log_x=self.params["scale_options.log_x"],
            log_y=self.params["scale_options.log_y"],
            log_z=self.params["scale_options.log_z"],
            color="color" if self.params["data.color.enabled"] else None,
            size="size" if self.params["data.size.enabled"] else None,
            labels=dict(
                x=self.params["data.x.label"],
                y=self.params["data.y.label"],
                z=self.params["data.z.label"],
                color=self.params["data.color.label"],
                size=self.params["data.size.label"]
            ),
            **colorscale
        )

        st.plotly_chart(fig)
