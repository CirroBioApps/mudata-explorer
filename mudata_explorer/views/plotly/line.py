import plotly.express as px
from streamlit.delta_generator import DeltaGenerator
from mudata_explorer.views.plotly.base import Plotly


class PlotlyLine(Plotly):

    type = "plotly-line"
    name = "Line Plot (Plotly)"
    help_text = "Display a series of data as a line graph using Plotly."
    schema = {
        "data": {
            "type": "dataframe",
            "columns": {
                "x": {"label": "x-axis"},
                "y": {"label": "y-axis"},
                "sort_by": {"label": "Sort By"},
                "color": {
                    "label": "Color",
                    "optional": True,
                    "colorscale": True
                },
            },
            "query": True,
        },
        "scale_options": {
            "type": "object",
            "label": "Scale Options",
            "properties": {
                "log_x": {
                    "type": "boolean",
                    "label": "Log Scale - X Axis"
                },
                "log_y": {
                    "type": "boolean",
                    "label": "Log Scale - Y Axis"
                }
            }
        }
    }

    def display(
        self,
        container: DeltaGenerator
    ):
        data, colorscale = self.fetch_dataframe("data")
        if data is None:
            container.write("Please select a data table")
            return
        data.sort_values("sort_by", inplace=True)

        fig = px.line(
            data,
            x="x",
            y="y",
            color="color" if self.params["data.color.enabled"] else None,
            log_x=self.params["scale_options.log_x"],
            log_y=self.params["scale_options.log_y"],
            labels=dict(
                x=self.params["data.x.label"],
                y=self.params["data.y.label"],
                color=self.params["data.color.label"]
            ),
            **colorscale
        )

        container.plotly_chart(fig)
