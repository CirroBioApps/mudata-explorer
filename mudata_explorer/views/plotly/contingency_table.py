from typing import Optional
import numpy as np
import pandas as pd
import plotly.express as px
from plotly.graph_objects import Figure
from plotly import graph_objects as go
import streamlit as st
from mudata_explorer.views.plotly.base import Plotly
from scipy.cluster import hierarchy
from scipy.stats import fisher_exact


class PlotlyContingencyTable(Plotly):

    type = "plotly-contingency-table"
    name = "Contingency Table (Plotly)"
    help_text = """
Compare the number of samples which are found in each combination of
two different columns of categories. This is a frequency table which
shows the number of times that each unique combination of values
is found in the data.

The primary purpose of this display is to identify when there is
a strong correlation between the values in two different columns.

The display can be used to show either:

- The number of items in each combination, or
- The odds ratio of each combination, calculated as the log2 of the
    ratio of the observed frequency to the expected frequency.

    """
    schema = {
        "data": {
            "type": "dataframe",
            "columns": {
                "x": {"label": "Category A"},
                "y": {"label": "Category B"}
            },
            "query": True,
            "sidebar": True
        },
        "formatting": {
            "type": "object",
            "label": "Formatting",
            "properties": {
                "colorscale": {
                    "type": "string",
                    "label": "Color Scale",
                    "help": "The color scale to use for the color metric.",
                    "default": "blues",
                    "enum": px.colors.named_colorscales()
                },
                "values": {
                    "type": "string",
                    "label": "Display Values",
                    "help": "The values to display in the table.",
                    "default": "Number of Items",
                    "enum": ["Number of Items", "Odds Ratio (log2)"],
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

    def build_figure(self) -> Figure:

        data: Optional[pd.DataFrame] = self.params.get("data.dataframe")
        if data is None:
            st.write("Please select two columns to compare")
            return

        # Make the contingency table
        table = pd.crosstab(
            data["y"],
            data["x"]
        )

        # Sort using linkage clustering
        table = self.sort_rows(
            self.sort_rows(table.T).T
        )

        # Calculate the enrichment
        prop = table / table.sum().sum()
        exp = pd.DataFrame({
            col: {
                row: prop.loc[row].sum() * prop[col].sum()
                for row in table.index
            }
            for col in table.columns
        })
        enrichment = pd.DataFrame({
            col: {
                row: (prop.loc[row, col] - exp.loc[row, col]) / exp.loc[row, col]
                for row in table.index
            }
            for col in table.columns
        }) * 100

        # Calculate the p-value for each cell
        pvals = pd.DataFrame({
            col: {
                row: fisher_exact(
                    [
                        [table.loc[row, col], table.loc[row].sum() - table.loc[row, col]],
                        [table.loc[:, col].sum() - table.loc[row, col], table.drop(columns=[col], index=[row]).sum().sum()]
                    ]
                )[1]
                for row in table.index
            }
            for col in table.columns
        })

        # Make the text for each cell
        disp = pd.DataFrame({
            col: {
                row: (
                    f"Observed: {table.loc[row, col]}<br>" +
                    f"Expected: {table.sum().sum() * exp.loc[row, col]:.1f}<br>" +
                    f"Enrichment: {enrichment.loc[row, col]:.2f}%<br>" +
                    f"p-value: {pvals.loc[row, col]:.2f}"
                ) if pvals.loc[row, col] < 0.1 else ""
                for row in table.index
            }
            for col in table.columns
        })

        # If the user has selected to display an odds ratio
        value_format = self.params["formatting.values"]
        if value_format == "Odds Ratio (log2)":
            table = pd.DataFrame({
                col: {
                    row: prop.loc[row, col] / exp.loc[row, col]
                    for row in table.index
                }
                for col in table.columns
            }).apply(np.log2)

        fig = go.Figure(
            data=go.Heatmap(
                z=table.values,
                x=table.columns,
                y=table.index,
                text=disp.values,
                texttemplate="%{text}",
                hovertext=table.applymap(
                    lambda i: (
                        f"{value_format}: {i}"
                    )
                ),
                zmid=0 if value_format == "Odds Ratio (log2)" else None,
                colorscale=self.params["formatting.colorscale"],
                colorbar_title=value_format
            )
        )
        fig.update_layout(
            xaxis_title=self.params["data.columns.x.label.value"],
            yaxis_title=self.params["data.columns.y.label.value"]
        )

        return fig
