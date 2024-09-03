import numpy as np
import pandas as pd
import plotly.express as px
import streamlit as st
from mudata_explorer.views.plotly.base import Plotly


class PlotlyScatterMulti(Plotly):

    type = "plotly-scatter-multi"
    name = "Scatterplot - Multiple Panels (Plotly)"
    help_text = """
Display a two dimensional distribution of data using Plotly.
Display the same set of points in multiple panels, each of which
is colored by a different column of data.

In addition to the data used for the x- and y-axes, the size
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
                "size": {"label": "size", "optional": True}
            },
            "query": True,
        },
        "colors": {
            "type": "dataframe",
            "label": "Color Columns",
            "help": "Each column will be used to color the points in a separate panel.", # noqa
            "select_columns": True,
            "query": "",
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
                "log_color": {
                    "type": "boolean",
                    "label": "Log Scale - Color Axis",
                    "sidebar": True
                }
            }
        },
        "formatting": {
            "type": "object",
            "label": "Formatting",
            "properties": {
                "colorscale": {
                    "type": "string",
                    "label": "Color Scale",
                    "help": "The color scale to use for the color metric.",
                    "default": "bluered",
                    "enum": px.colors.named_colorscales()
                },
                "opacity": {
                    "type": "float",
                    "label": "Point Opacity",
                    "default": 1.0,
                    "min_value": 0.,
                    "max_value": 1.,
                    "step": 0.1,
                    "sidebar": True
                },
                "ncols": {
                    "type": "integer",
                    "label": "Number of Columns",
                    "default": 1,
                    "min_value": 1
                },
                "color_label": {
                    "type": "string",
                    "label": "Color Label",
                    "default": "Abundance"
                }
            }
        }
    }

    def display(self):

        # Get the tables provided by the user
        data = self.params.get("data.dataframe")
        if data is None:
            st.write("Please select a data table")
            return
        colors = self.params.get("colors.dataframe")
        if colors is None:
            st.write("Please select a colors table")
            return

        # Make sure that the indices are the same
        index = data.index.intersection(colors.index)
        assert len(index) > 0, \
            "No common indices found between the data and colors tables."
        data = data.loc[index]
        colors = colors.loc[index]

        # Merge the tables by making new columns 'color_val'
        # which has the values from the color columns
        # as well as "color_label" which has the column name.
        # The resulting DataFrame will be long, with a new
        # copy for each distinct column
        data_long = pd.concat([
            data.assign(
                color_val=color_val,
                color_label=color_label
            )
            for color_label, color_val in colors.items()
        ])

        # If the log scale has been set for the colors
        if self.params["scale_options.log_color"]:
            data_long = (
                data_long
                .query("color_val > 0")
                .assign(
                    color_val=lambda df: df["color_val"].apply(np.log10)
                )
            )

        try:
            opacity = float(self.params["formatting.opacity"])
        except ValueError:
            opacity = 1.0

        fig = px.scatter(
            data_long.reset_index(),
            x="x",
            y="y",
            hover_name=(
                "index" if data.index.name is None else
                data.index.name
            ),
            log_x=self.params["scale_options.log_x"],
            log_y=self.params["scale_options.log_y"],
            color="color_val",
            facet_col="color_label",
            facet_col_wrap=int(self.params["formatting.ncols"]),
            size="size" if self.params["data.columns.size.enabled.value"] else None,
            opacity=opacity,
            labels=dict(
                x=self.params["data.columns.x.label.value"],
                y=self.params["data.columns.y.label.value"],
                size=self.params["data.columns.size.label.value"],
                color_val=self.params["formatting.color_label"]
            ),
            color_continuous_scale=self.params["formatting.colorscale"]
        )
        fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))

        st.plotly_chart(fig)
