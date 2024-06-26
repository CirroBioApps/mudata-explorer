import pandas as pd
import plotly.express as px
from streamlit.delta_generator import DeltaGenerator
from mudata_explorer.views.plotly.base import Plotly


class PlotlyBoxMulti(Plotly):

    type = "plotly-box-multiple"
    name = "Box Plot - Multiple Measurements (Plotly)"
    help_text = """
    Display multiple columns of data as a box graph using Plotly, summarizing
    the data in terms of the median, quartiles, and outliers.

    A collection of columns are used to define the values on the y-axis, and a
    second column is used for the categorical groups which are
    displayed on the x-axis.
    """
    schema = {
        "table": {
            "type": "object",
            "label": "Data Table",
            "properties": {
                "data": {
                    "type": "dataframe",
                    "label": "Data",
                    "help": "Select the measurement values to summarize",
                    "select_columns": True,
                    "query": "",
                },
                "category": {
                    "type": "dataframe",
                    "label": "Category",
                    "help": "Select the column containing category labels", # noqa
                    "columns": {"category": {"label": "Category"}}
                }
            }
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
        "display_options": {
            "type": "object",
            "label": "Display Options",
            "properties": {
                "ncols": {
                    "type": "integer",
                    "label": "Number of Columns",
                    "default": 1,
                    "min_value": 1
                },
                "outliers": {
                    "type": "boolean",
                    "label": "Show Outliers",
                    "default": True
                }
            }
        }
    }

    def display(self, container: DeltaGenerator):

        data: pd.DataFrame = self.params.get("table.data.dataframe")
        if data is None:
            container.write("Please select a data table")
            return
        category: pd.Series = (
            self.params
            ["table.category.dataframe"]
            ["category"]
        )

        # Get the shared indices
        index = data.index.intersection(category.index)

        # Make sure that there is some degree of intersection
        msg = "No common indices found between the data and category tables."
        assert len(index) > 0, msg

        # Subset each to just those indices
        data = data.loc[index]
        category = category.loc[index]
        data.index.name = "index"

        # Make a long DataFrame which has all of the values
        data_long = (
            data
            .rename_axis(columns=None)
            .reset_index()
            .melt(id_vars="index")
            .assign(
                category=lambda df: df["index"].apply(category.get)
            )
        )

        fig = px.box(
            data_long,
            x="category",
            y="value",
            facet_col=(
                "variable"
                if "variable" in data_long.columns
                else "var"
            ),
            boxmode="overlay",
            facet_col_wrap=int(self.params["display_options.ncols"]),
            log_y=self.params["scale_options.log_y"],
            labels=dict(
                category=self.params["table.category.category.label"]
            ),
            points=(
                "outliers"
                if self.params["display_options.outliers"]
                else False
            )
        )
        fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
        fig.update_yaxes(matches=None)

        container.plotly_chart(fig)
