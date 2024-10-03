import plotly.express as px
import plotly.graph_objects as go
import streamlit as st
from mudata_explorer.views.plotly.base import Plotly
from scipy.stats import linregress


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
        "stats_options": {
            "type": "object",
            "label": "Statistics Options",
            "properties": {
                "linregress": {
                    "type": "boolean",
                    "label": "Show Linear Regression (scipy.stats.linregress)",
                    "default": False,
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
                },
                "legend": {
                    "type": "string",
                    "label": "Legend",
                    "multiline": True
                }
            }
        }
    }

    def build_figure(self) -> go.Figure:

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

        # Set the default title
        title = self.params.get("formatting.title", "")

        # If the user has enabled linear regression
        if self.params["stats_options.linregress"]:

            # Run the linear regression
            slope, intercept, r_value, p_value, std_err = linregress(
                data["x"],
                data["y"]
            )

            # Add the r_value and p_value to the title
            if title:
                title += f" (Linear Regression: r={r_value:.2f}, p={p_value:.2e})"
            else:
                title = f"Linear Regression: r={r_value:.2f}, p={p_value:.2e}"

        # Make the figure
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
            title=title,
            **colorscale
        )

        # If the user has enabled linear regression
        if self.params["stats_options.linregress"]:

            # Add the linear regression line to the figure behind other elements
            fig.add_trace(
                go.Scatter(
                    x=[data["x"].min(), data["x"].max()],
                    y=[slope * data["x"].min() + intercept, slope * data["x"].max() + intercept],
                    mode="lines",
                    line=dict(color="grey", dash="dash"),
                    name="Linear Regression"
                )
            )
            # Move the linear regression line trace to the first position
            fig.data = fig.data[-1:] + fig.data[:-1]

        return fig
