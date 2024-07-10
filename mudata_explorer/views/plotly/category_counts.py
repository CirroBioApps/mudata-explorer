import plotly.express as px
from streamlit.delta_generator import DeltaGenerator
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
                    "label": "Chi-Square Contingency Test"
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
            container.error("Color scale must be categorical for this display")

        # Count up the number of rows for each unique value
        df = (
            data
            .groupby(
                ["x", "color"]
                if self.params["data.color.enabled"] else "x"
            )
            .apply(len)
            .rename("count")
            .reset_index()
        )

        if self.params["annotation_options.chisquare"]:
            if self.params["data.color.enabled"]:
                chi2, p, _, _ = chi2_contingency(
                    df.pivot(index="x", columns="color", values="count").fillna(0)
                )
                if p > 0.001:
                    colorscale['title'] = f"Chi-Square: {chi2:.2f}, p-value: {p:.4f}"
                else:
                    colorscale['title'] = f"Chi-Square: {chi2:.2f}, p-value: {p:.2E}"

        fig = px.bar(
            df,
            x="x",
            y="count",
            log_y=self.params["scale_options.log_y"],
            color="color" if self.params["data.color.enabled"] else None,
            labels=dict(
                x=self.params["data.x.label"],
                count=self.params["ylabel"],
                color=self.params["data.color.label"]
            ),
            barmode=self.params["barmode"],
            text_auto=self.params.get("annotation_options.show_values", False) is True,
            **colorscale
        )

        container.plotly_chart(fig)
