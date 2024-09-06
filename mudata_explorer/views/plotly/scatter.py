import plotly.express as px
import streamlit as st
from mudata_explorer.views.plotly.base import Plotly


class PlotlyScatter(Plotly):

    type = "plotly-scatter"
    name = "Scatterplot (Plotly)"
    help_text = """
Display a two dimensional distribution of data using Plotly.

In addition to the data used for the x- and y-axes, the size and color
of the points can be set to represent additional dimensions of the data.

Additional formatting options include setting the opacity of the points
and using a log scale for the x- and y-axes.
"""
    schema = {
        "data": {
            "type": "dataframe",
            "columns": {
                "x": {"label": "x-axis"},
                "y": {"label": "y-axis"},
                "size": {"label": "size", "optional": True},
                "color": {
                    "label": "Color",
                    "optional": True,
                    "colorscale": True
                },
            },
            "query": True,
            "sidebar": True
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
                }
            }
        },
        "formatting": {
            "type": "object",
            "label": "Formatting",
            "properties": {
                "opacity": {
                    "type": "float",
                    "label": "Point Opacity",
                    "default": 1.0,
                    "min_value": 0.,
                    "max_value": 1.,
                    "step": 0.1,
                    "sidebar": True
                },
                "title": {
                    "type": "string",
                    "label": "Figure Title",
                    "default": "",
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
        if "x" not in data.columns or "y" not in data.columns:
            st.write("Please select columns for the x and y axes")
            return
        try:
            opacity = float(self.params["formatting.opacity"])
        except ValueError:
            opacity = 1.0

        fig = px.scatter(
            data.reset_index(),
            x="x",
            y="y",
            hover_name=(
                "index" if data.index.name is None else
                data.index.name
            ),
            log_x=self.params["scale_options.log_x"],
            log_y=self.params["scale_options.log_y"],
            color="color" if self.params["data.columns.color.enabled.value"] else None,
            size="size" if self.params["data.columns.size.enabled.value"] else None,
            opacity=opacity,
            labels=dict(
                x=self.params["data.columns.x.label.value"],
                y=self.params["data.columns.y.label.value"],
                color=self.params["data.columns.color.label.value"],
                size=self.params["data.columns.size.label.value"]
            ),
            title=self.params["formatting.title"],
            **colorscale
        )

        st.plotly_chart(fig)
