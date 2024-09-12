import pandas as pd
import plotly.express as px
import streamlit as st
from mudata_explorer.views.plotly.base import Plotly
from scipy.cluster import hierarchy
from typing import Optional


class PlotlyCategorySummarizeValues(Plotly):

    type = "plotly-category-summarize-values"
    name = "Category Summary (Plotly)"
    help_text = """
    Summarize the values in a collection of columns,
    grouped by a single column which contains categories.

    The plot will be laid out with the categories on the y-axis,
    and each of the columns will be shown on the x-axis.

    For each of the categories, the values in each of the columns
    will be summarized in terms of:

    - The mean value
    - The median value
    - The number of positive, non-zero values
    - The number of non-null values

    The size and color of the points can be set to represent
    any of those summary statistics.
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
        "formatting": {
            "type": "object",
            "label": "Formatting",
            "properties": {
                "size": {
                    "type": "string",
                    "enum": ["Mean", "Median", "Positive", "Non-Null", "None"],
                    "default": "Mean",
                    "label": "Point Size",
                    "help": "Metric used to scale the size of each point.",
                    "sidebar": True
                },
                "color": {
                    "type": "string",
                    "enum": ["Mean", "Median", "Positive", "Non-Null", "None"],
                    "default": "Non-Null",
                    "label": "Point Color",
                    "help": "Metric used to scale the color of each point."
                },
                "colorscale": {
                    "type": "string",
                    "label": "Color Scale",
                    "help": "The color scale to use for the color metric.",
                    "default": "bluered",
                    "enum": px.colors.named_colorscales()
                },
                "sort_by": {
                    "type": "string",
                    "label": "Sort By",
                    "enum": ["Labels", "Values"],
                    "default": "Labels",
                    "sidebar": True
                },
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
                "transpose": {
                    "type": "boolean",
                    "label": "Transpose",
                    "default": False,
                    "help": "Swap the x and y axes."
                }
            }
        }
    }

    @staticmethod
    def sort_rows(df: pd.DataFrame) -> pd.DataFrame:
        if df.shape[0] <= 1:
            return df
        return df.reindex(
            index=df.index[
                hierarchy.leaves_list(
                    hierarchy.linkage(
                        df.values,
                        method="ward"
                    )
                )
            ]
        )

    def display(self):

        data: Optional[pd.DataFrame] = self.params.get("table.data.dataframe")
        if data is None:
            st.write("Please select a data table")
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

        # Calculate the summary metrics
        summary = (
            data
            .groupby(category)
            .apply(
                lambda group_df: group_df.apply(
                    lambda x: pd.Series({
                        "Mean": x.mean(),
                        "Median": x.median(),
                        "Positive": (x.dropna() > 0).sum() / x.shape[0],
                        "Non-Null": x.notnull().sum() / x.shape[0]
                    })
                )
            )
            .reset_index()
            .rename(columns={"level_1": "metric"})
            .melt(id_vars=["category", "metric"], var_name="column")
            .pivot_table(
                index=["column", "category"],
                columns="metric",
                values="value"
            )
        )

        # Add scaled values for the size and color
        summary = pd.concat(
            [
                summary,
                summary.apply(
                    lambda c: (
                        (c - c.min()) / (c.max() - c.min()) + 0.01
                        if c.min() < c.max()
                        else c
                    )
                ).rename(
                    columns=lambda c: f"{c}-scaled"
                )
            ],
            axis=1
        )

        labels = dict()
        if self.params['formatting.size'] != "None":
            size = self.params['formatting.size'] + "-scaled"
            labels[size] = self.params['formatting.size']
        else:
            size = None

        if self.params['formatting.color'] != "None":
            color = self.params['formatting.color'] + "-scaled"
            labels[color] = self.params['formatting.color']
        else:
            color = None

        summary.reset_index(inplace=True)

        # If all of the columns start with the same prefix
        if (
            summary["column"].apply(lambda x: ":" in x).all()
            and summary["column"].apply(lambda x: x.split(":")[0]).nunique() == 1
        ):
            summary["column"] = summary["column"].apply(
                lambda x: x.split(":", 1)[1]
            )

        if self.params['formatting.sort_by'] == "Labels":
            category_orders = {
                cname: summary[cname].drop_duplicates().tolist()
                for cname in ["column", "category"]
            }
        else:
            wide_df = summary.pivot_table(
                index="category",
                columns="column",
                values=size
            )
            category_orders = dict(
                category=self.sort_rows(wide_df).index,
                column=self.sort_rows(wide_df.T).index
            )

        if self.params["formatting.transpose"]:
            x = "category"
            y = "column"
        else:
            x = "column"
            y = "category"

        fig = px.scatter(
            summary,
            x=x,
            y=y,
            size=size,
            color=color,
            color_continuous_scale=self.params['formatting.colorscale'],
            hover_data=["Mean", "Median", "Positive", "Non-Null"],
            category_orders=category_orders,
            labels=dict(
                column=(
                    "Observation"
                    if self.params["table.data.axis"]
                    else "Variable"
                ),
                category=self.params["table.category.columns.category.label"],
                **labels
            ),
            title=self.params["formatting.title"]
        )

        st.plotly_chart(fig)
        self.show_legend(key="formatting.legend")
