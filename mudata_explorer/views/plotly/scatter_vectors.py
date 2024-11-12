import plotly.express as px
from plotly.graph_objects import Figure
import pandas as pd
import streamlit as st
from mudata_explorer.views.plotly.base import Plotly


class PlotlyScatterVectors(Plotly):

    type = "plotly-scatter-vectors"
    name = "Scatterplot with Vectors (Plotly)"
    help_text = """
Display a two dimensional distribution of data as a scatter plot,
while also overlaying a set of vectors on the plot.

The length of the vectors will be scaled to the span of the points.
The name of each vector may optionally be displayed next to the vector.

In addition to the data used for the x- and y-axes, the size and color
of the points can be set to represent additional dimensions of the data.

Additional formatting options include setting the opacity of the points
and using a log scale for the x- and y-axes.
"""
    schema = {
        "data": {
            "type": "dataframe",
            "label": "Data used to place the points on the plot",
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
            "sidebar": False
        },
        "vectors": {
            "type": "dataframe",
            "label": "Data used to place the vectors on the plot",
            "columns": {
                "x": {"label": "x-axis"},
                "y": {"label": "y-axis"},
                "label": {
                    "label": "Vector Label",
                    "optional": True
                }
            },
            "query": True,
            "sidebar": False
        },
        "scale_options": {
            "type": "object",
            "label": "Scale Options",
            "properties": {
                "log_x": {
                    "type": "boolean",
                    "label": "Log Scale - X Axis",
                    "sidebar": False
                },
                "log_y": {
                    "type": "boolean",
                    "label": "Log Scale - Y Axis",
                    "sidebar": False
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
                "vector_opacity": {
                    "type": "float",
                    "label": "Vector Opacity",
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

    def build_figure(self) -> Figure:

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

        vectors, _ = self.fetch_dataframe("vectors")
        if vectors is None:
            st.write("Please select a data table for the vectors")
            return
        if "x" not in vectors.columns or "y" not in vectors.columns:
            st.write("Please select columns for the x and y axes of the vectors")
            return
        
        if "color_discrete_sequence" not in colorscale:
            colorscale["color_discrete_sequence"] = px.colors.qualitative.Plotly

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

        # Scale the range of the vectors to the range of the points
        vectors = _scale_vectors(data, vectors)

        try:
            vector_opacity = float(self.params["formatting.vector_opacity"])
        except ValueError:
            vector_opacity = 1.0

        for i, row in vectors.iterrows():
            # Add the arrow
            fig.add_annotation(
                x=row["x"], xref="x", yref="y",
                y=row["y"], axref="x", ayref="y",
                ax=0,
                ay=0,
                showarrow=True,
                opacity=vector_opacity
            )
            # Add the label
            fig.add_annotation(
                x=row["x"], xref="x", yref="y",
                y=row["y"], axref="x", ayref="y",
                text=row["label"] if "label" in vectors.columns else i,
                showarrow=False,
                opacity=vector_opacity
            )

        return fig


def _scale_vectors(data: pd.DataFrame, vectors: pd.DataFrame):
    """Scale the vectors to the range of the data."""

    # For each dimension
    for kw in ["x", "y"]:
        # First scale to 0-1
        vectors[kw] = (vectors[kw] - vectors[kw].min()) / _range(vectors[kw])

        # Now align to the range of the data
        vectors[kw] = vectors[kw] * _range(data[kw]) + data[kw].min()

    return vectors
    

def _range(vec: pd.Series):
    return vec.max() - vec.min()
