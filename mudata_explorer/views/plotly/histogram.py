import plotly.express as px
from plotly.graph_objects import Figure
from scipy.stats import f_oneway, kruskal
import streamlit as st
from mudata_explorer.views.plotly.base import Plotly


class PlotlyHistogram(Plotly):

    type = "plotly-histogram"
    name = "Histogram (Plotly)"
    help_text = """
Display a distribution of values as a frequency histogram using Plotly.

Optionally include a grouping column to display multiple histograms
which are overlaid on the same plot using different colors.
    """
    schema = {
        "data": {
            "type": "dataframe",
            "columns": {
                "value": {"label": "Value"},
                "grouping": {
                    "label": "Grouping",
                    "optional": True,
                    "colorscale": True,
                    "is_categorical": True
                }
            },
            "query": True
        },
        "scale_options": {
            "type": "object",
            "label": "Scale Options",
            "properties": {
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
                "title": {
                    "type": "string",
                    "label": "Title",
                    "default": "",
                    "sidebar": True
                },
                "legend": {
                    "type": "string",
                    "label": "Legend",
                    "multiline": True
                },
                "nbins": {
                    "type": "integer",
                    "label": "Number of Bins",
                    "default": 20,
                    "min_value": 1,
                    "sidebar": True
                }
            }
        },
        "statistics": {
            "type": "object",
            "label": "Statistics",
            "properties": {
                "compare_groups": {
                    "type": "string",
                    "label": "Compare Values Between Groups",
                    "default": "Disabled",
                    "enum": ["Disabled", "ANOVA", "Kruskal-Wallis"],
                    "sidebar": True
                }
            }
        }
    }

    def build_figure(self) -> Figure:

        data, colorscale = self.fetch_dataframe("data")
        if data is None:
            st.write("Please select a data table")
            return
        if "value" not in data.columns:
            st.write("Please select a column to use as the value")
            return

        if "grouping" in data.columns and "color_continuous_scale" in colorscale:
            st.error("Color scale must be categorical for box plots.")
            return

        title = self.params["formatting.title"]
        nbins = self.params["formatting.nbins"]
        log_y = self.params["scale_options.log_y"]

        method = self.params["statistics.compare_groups"]
        if method != "Disabled" and "grouping" in data.columns:
            vals = [
                d["value"].dropna().values for group, d in data.groupby("grouping")
                if d["value"].notnull().sum() > 1 and group is not None
            ]
            title = self._add_stats_title(title, method, vals)

        if "color_discrete_sequence" not in colorscale:
            colorscale["color_discrete_sequence"] = px.colors.qualitative.Plotly

        fig = px.histogram(
            data,
            x="value",
            nbins=nbins,
            log_y=log_y,
            color="grouping" if "grouping" in data.columns else None,
            labels=dict(
                value=self.params["data.columns.value.label.value"],
                grouping=self.params["data.columns.grouping.label.value"]
            ),
            title=title,
            **colorscale
        )

        return fig
