import pandas as pd
import plotly.express as px
from streamlit.delta_generator import DeltaGenerator
from mudata_explorer.views.plotly.base import Plotly
from scipy.cluster import hierarchy


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
                    "help": "Metric used to scale the size of each point."
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
                }
            }
        }
    }

    @staticmethod
    def sort_rows(df: pd.DataFrame) -> pd.DataFrame:
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

        if self.params['formatting.sort_by'] == "Labels":
            category_orders = {
                cname: summary.reset_index()[cname].drop_duplicates().tolist()
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

        fig = px.scatter(
            summary.reset_index(),
            x="column",
            y="category",
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
                category=self.params["table.category.category.label"],
                **labels
            )
        )

        container.plotly_chart(fig)
