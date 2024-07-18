import plotly.express as px
from scipy.stats import f_oneway, kruskal
from streamlit.delta_generator import DeltaGenerator
from mudata_explorer.views.plotly.base import Plotly


class PlotlyBox(Plotly):

    type = "plotly-box"
    name = "Box Plot - Single Measurement (Plotly)"
    help_text = """
    Display a single column of data as a box graph using Plotly, summarizing
    the data in terms of the median, quartiles, and outliers.

    A single column is used to define the values on the y-axis, and a
    second column is used for the categorical groups which are
    displayed on the x-axis.
    """
    schema = {
        "data": {
            "type": "dataframe",
            "columns": {
                "x": {"label": "x-axis"},
                "y": {"label": "y-axis"},
                "color": {
                    "label": "Color",
                    "optional": True,
                    "colorscale": True,
                    "is_categorical": True
                },
            },
            "query": True,
        },
        "scale_options": {
            "type": "object",
            "label": "Scale Options",
            "properties": {
                "log_y": {
                    "type": "boolean",
                    "label": "Log Scale - Y Axis"
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
                    "default": ""
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
                    "enum": ["Disabled", "ANOVA", "Kruskal-Wallis"]
                }
            }
        }
    }

    def display(self, container: DeltaGenerator):

        data, colorscale = self.fetch_dataframe("data")
        if data is None:
            container.write("Please select a data table")
            return

        if "color_continuous_scale" in colorscale:
            container.error("Color scale must be categorical for box plots.")

        title = self.params["formatting.title"]

        method = self.params["statistics.compare_groups"]
        if method != "Disabled":
            vals = [d["y"].values for _, d in data.groupby("x")]
            if method == "ANOVA":
                res = f_oneway(*vals)
            else:
                assert method == "Kruskal-Wallis"
                res = kruskal(*vals)
            # Format the p-value as a string
            pvalue = (
                f"{res.pvalue:.4f}"
                if res.pvalue > 0.0001
                else f"{res.pvalue:.2e}"
            )

            # Add it to the title
            formatted_result = f"{method} p-value: {pvalue}"
            if len(title) > 0:
                title = f"{title} ({formatted_result})"
            else:
                title = formatted_result

        fig = px.box(
            data,
            x="x",
            y="y",
            log_y=self.params["scale_options.log_y"],
            color="color" if self.params["data.color.enabled"] else None,
            labels=dict(
                x=self.params["data.x.label"],
                y=self.params["data.y.label"],
                color=self.params["data.color.label"]
            ),
            title=title,
            **colorscale
        )

        container.plotly_chart(fig)
