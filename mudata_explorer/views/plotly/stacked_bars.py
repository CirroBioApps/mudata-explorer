import pandas as pd
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import streamlit as st
from mudata_explorer.views.plotly.base import Plotly
from scipy.cluster import hierarchy
from typing import Optional


class PlotlyStackedBars(Plotly):

    type = "plotly-stacked-bars"
    name = "Stacked Bar Chart (Plotly)"
    help_text = """
Show the values in a collection of columns, optionally
annotated by a single column which contains categories.
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
                    "columns": {
                        "category": {
                            "label": "Category"
                        }
                    },
                    "optional": True
                }
            }
        },
        "formatting": {
            "type": "object",
            "label": "Formatting",
            "properties": {
                "max_features": {
                    "type": "integer",
                    "label": "Maximum Number of Features",
                    "optional": True
                },
                "sort_cols_by": {
                    "type": "string",
                    "label": "Sort Features By",
                    "enum": ["Labels", "Values"],
                    "default": "Labels",
                    "sidebar": True
                },
                "sort_rows_by": {
                    "type": "string",
                    "label": "Sort Bars By",
                    "enum": ["Labels", "Values", "Category"],
                    "default": "Labels",
                    "sidebar": True
                },
                "title": {
                    "type": "string",
                    "label": "Title",
                    "default": "",
                    "sidebar": True
                },
                "yaxis_title": {
                    "type": "string",
                    "label": "Y Axis Title",
                    "default": ""
                },
                "feature_label": {
                    "type": "string",
                    "label": "Feature",
                    "default": ""
                },
                "max_y": {
                    "type": "float",
                    "label": "Maximum Y Value",
                    "help": "If specified, set the maximum value of the y axis",
                    "default": 0
                }
            }
        }
    }
    _below_threshold_label = "(Below Threshold)"

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

        # Group together all columns beyond the `max_features` threshold
        max_features = self.params.get("formatting.max_features.value")
        if max_features is not None and max_features > 0 and data.shape[1] > max_features:
            feature_rank = data.mean().sort_values(ascending=False)
            to_drop = feature_rank.iloc[max_features:].index
            data = (
                data
                .drop(columns=to_drop)
                .assign(**{
                    self._below_threshold_label: data.loc[:, to_drop].sum(axis=1)
                })
            )

        # If no category was given
        if self.params.get("table.category.dataframe") is None:
            category = None
        # If a category DataFrame was provided
        else:
            category_df: pd.DataFrame = self.params["table.category.dataframe"]

            # If a column was provided
            if "category" in category_df.columns.values:
                category: pd.Series = category_df["category"]
            else:
                category = None

        if category is not None:

            # Get the shared indices
            index = data.index.intersection(category.index)

            # Make sure that there is some degree of intersection
            msg = "No common indices found between the data and category tables."
            assert len(index) > 0, msg

            # Subset each to just those indices
            data = data.loc[index]
            category = category.loc[index]

        # Sort the table
        if self.params['formatting.sort_rows_by.value'] == "Labels":
            data = data.sort_index()
        elif self.params['formatting.sort_rows_by.value'] == "Values":
            data = self.sort_rows(data)
        else:
            assert self.params['formatting.sort_rows_by.value'] == "Category"
            assert category is not None, "Must provide a category to sort by"
            data = pd.concat([
                self.sort_rows(d)
                for _, d in data.groupby(category)
            ])

        if self.params['formatting.sort_cols_by.value'] == "Labels":
            data = data.sort_index(axis=1)
        else:
            assert self.params['formatting.sort_cols_by.value'] == "Values"
            data = self.sort_rows(data.T).T

        # Set up the layout arguments
        yaxis_title = self.params.get("formatting.yaxis_title.value")
        feature_label = self.params.get("formatting.feature_label.value")
        layout_args = dict(
            barmode='stack',
            title=self.params.get("formatting.title.value"),
            yaxis_title=yaxis_title
        )

        # If the maximum Y value was specified
        max_y = self.params.get("formatting.max_y.value", 0)
        if max_y > 0:
            layout_args["yaxis_range"] = [0, max_y]

        if category is None:
            # Just make the stacked bars
            fig = go.Figure(data=self._make_bars(data, feature_label))
        else:
            # Make the stacked bars with categories
            fig = make_subplots(
                rows=2,
                cols=1,
                shared_xaxes=True,
                vertical_spacing=0.01,
                row_heights=[0.9, 0.1]
            )
            fig.add_traces(self._make_bars(data, feature_label), rows=1, cols=1)
            category_name = self.params.get("table.category.category.label.value", "")
            fig.add_traces(self._make_annot(category, category_name), rows=2, cols=1)
            # Set the ticks on yaxis2 to be invisible
            fig.update_yaxes(
                tickvals=[0.5],
                ticktext=[category_name],
                row=2,
                col=1
            )

        # Change the bar mode and add titles
        fig.update_layout(**layout_args)

        st.plotly_chart(fig)

    def _make_bars(self, data: pd.DataFrame, yaxis_title: str):
        return [
            go.Bar(
                name=str(col),
                x=data.index,
                y=data[col],
                legendgroup=yaxis_title,
                legendgrouptitle_text=yaxis_title
            )
            for col in (
                [self._below_threshold_label]
                if self._below_threshold_label in data.columns
                else []
            ) + [
                cname
                for cname in data.columns[::-1]
                if cname != self._below_threshold_label
            ]
        ]

    def _make_annot(self, category: pd.Series, category_name: str):
        return [
            go.Bar(
                name=str(name),
                x=category.index[category == name],
                y=[1 for _ in range((category == name).sum())],
                hovertext=str(name),
                legendgroup=category_name,
                legendgrouptitle_text=category_name
            )
            for name in category.unique()
        ]
