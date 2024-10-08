import plotly.express as px
from plotly.graph_objects import Figure
import streamlit as st
from mudata_explorer.views.plotly.base import Plotly
from scipy.stats import chi2_contingency


class PlotlyCategoryCount(Plotly):

    type = "plotly-category-count"
    name = "Category Count (Plotly)"
    help_text = """
    Show a bar graph which depicts the values in a single column,
    showing the number of times that each unique value is present.
    This is a frequency histogram showing categorical values.

    Optionally, a second column can be used to color the bars.
    When the color column is used, the position may be set to either
    stack vertically, or to group the bars side by side.
    """
    schema = {
        "data": {
            "type": "dataframe",
            "columns": {
                "x": {"label": "x-axis"},
                "color": {
                    "label": "Color",
                    "optional": True,
                    "colorscale": True,
                    "is_categorical": True
                },
            },
            "query": True,
        },
        "ylabel": {
            "type": "string",
            "label": "Y Axis Label",
            "help_text": "The label to use for the y-axis.",
            "default": "Count"
        },
        "barmode": {
            "type": "string",
            "label": "Bar Color Mode",
            "help_text": "If a color column is used, should the bars be stacked or grouped?", # noqa
            "enum": ["stack", "group"],
            "sidebar": True
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
                }
            }
        },
        "annotation_options": {
            "type": "object",
            "label": "Annotation Options",
            "properties": {
                "show_values": {
                    "type": "boolean",
                    "label": "Show Values"
                },
                "chisquare": {
                    "type": "boolean",
                    "label": "Chi-Square Contingency Test",
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

        # If a continuous color scale was selected
        if "color_continuous_scale" in colorscale:
            # Default to a categorical color scale
            colorscale = dict(color_discrete_sequence=px.colors.qualitative.D3)

        # Count up the number of rows for each unique value
        df = (
            data
            .groupby(
                ["x", "color"]
                if self.params["data.columns.color.enabled.value"] else "x"
            )
            .apply(len)
            .rename("count")
            .reset_index()
        )

        title = self.params["formatting.title"]

        if self.params["annotation_options.chisquare"]:
            if self.params["data.columns.color.enabled.value"]:
                chi2, p, _, _ = chi2_contingency(
                    df.pivot(index="x", columns="color", values="count").fillna(0)
                )
                if p > 0.001:
                    chi2_res = f"Chi-Square: {chi2:.2f}, p-value: {p:.4f}"
                else:
                    chi2_res = f"Chi-Square: {chi2:.2f}, p-value: {p:.2E}"

                if len(title) > 0:
                    title = f"{title} ({chi2_res})"
                else:
                    title = chi2_res

        fig = px.bar(
            df,
            x="x",
            y="count",
            log_y=self.params["scale_options.log_y"],
            color="color" if self.params["data.columns.color.enabled.value"] else None,
            labels=dict(
                x=self.params["data.columns.x.label.value"],
                count=self.params["ylabel.value"],
                color=self.params["data.columns.color.label.value"]
            ),
            barmode=self.params["barmode.value"],
            text_auto=self.params.get("annotation_options.show_values.value", False) is True,
            title=title,
            **colorscale
        )

        return fig
