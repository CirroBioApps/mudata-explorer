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
                    "columns": {"category": {"label": "Category"}},
                    "optional": True,
                    "query": ""
                }
            }
        },
        "variable_options": {
            "type": "object",
            "label": "Variable Options",
            "properties": {
                "axis": {
                    "type": "string",
                    "label": "Show Value On",
                    "enum": ["X-Axis", "Y-Axis"],
                    "default": "Y-Axis",
                },
                "log_values": {
                    "type": "boolean",
                    "label": "Log Scale"
                },
                "sort_by": {
                    "type": "string",
                    "enum": ["Mean", "Median", "Name"],
                    "default": "Mean",
                    "label": "Sort By",
                    "help": "How to sort the variables"
                }
            }
        },
        "category_options": {
            "type": "object",
            "label": "Category Options",
            "properties": {
                "axis": {
                    "type": "string",
                    "enum": ["Axis", "Facet", "Color"],
                    "default": "Axis",
                    "label": "Show Category As",
                    "help": "How to display the category data"
                },
                "sort_by": {
                    "type": "string",
                    "enum": ["Mean", "Median", "Name"],
                    "default": "Mean",
                    "label": "Sort By",
                    "help": "How to sort the categories"
                }
            }
        },
        "display_options": {
            "type": "object",
            "label": "Display Options",
            "properties": {
                "ncols": {
                    "type": "integer",
                    "label": "Number of Columns (optional)",
                    "default": 1,
                    "min_value": 1
                },
                "outliers": {
                    "type": "object",
                    "label": "Show Outlier Points",
                    "properties": {
                        "enabled": {
                            "type": "boolean",
                            "label": "Enabled",
                            "default": True
                        }
                    }
                },
                "title": {
                    "type": "string",
                    "label": "Title",
                    "default": ""
                },
                "var_label": {
                    "type": "string",
                    "label": "Variable Label",
                    "default": "Variable"
                },
                "val_label": {
                    "type": "string",
                    "label": "Value Label",
                    "default": "Value"
                },
                "height": {
                    "type": "integer",
                    "label": "Height",
                    "default": 500,
                    "min_value": 100
                }
            }
        }
    }

    def display(self, container: DeltaGenerator):

        data: pd.DataFrame = self.params.get("table.data.dataframe")
        if data is None:
            container.write("Please select a data table")
            return

        category = self.params.get("table.category.dataframe")
        if category is not None:
            category: pd.Series = category["category"]

            # Get the shared indices
            index = data.index.intersection(category.index)

            # Make sure that there is some degree of intersection
            msg = "No common indices found between the data and category tables."
            assert len(index) > 0, msg

            # Subset each to just those indices
            data = data.loc[index]
            category = category.loc[index]

        # Note that at this point the category may be None, in which case
        # we will format the plot differently

        # Clean up the index name
        data.index.name = "index"

        # Make a long DataFrame which has all of the values
        data_long = (
            data
            .rename_axis(columns=None)
            .reset_index()
            .melt(id_vars="index")
        )

        # Set up the elements which are always the same
        kwargs = dict(
            y="value",
            facet_col_wrap=int(self.params["display_options.ncols"]),
            height=self.params["display_options.height"],
            title=self.params["display_options.title"]
        )
        var_name = (
            "variable"
            if "variable" in data_long.columns
            else "var"
        )

        # If the category was provided
        if category is not None:

            # Add it to the table
            data_long = data_long.assign(
                category=lambda df: df["index"].apply(category.get)
            )

            # Format the display with the category as a facet
            if self.params["category_options.axis"] == "Facet":
                kwargs['x'] = var_name
                kwargs['facet_col'] = "category"
                kwargs['boxmode'] = "overlay"

            # Show the category on the x-axis
            elif self.params["category_options.axis"] == "Axis":
                kwargs['x'] = "category"
                kwargs['color'] = var_name
            else:
                # Display as a color
                assert self.params["category_options.axis"] == "Color"
                kwargs['x'] = var_name
                kwargs['color'] = "category"

        # If no category was provided
        else:
            # Format the display with no category
            kwargs["x"] = var_name

        # If the user wants to show the variable on the x-axis
        if self.params["variable_options.axis"] == "X-Axis":
            kwargs["x"], kwargs["y"] = kwargs["y"], kwargs["x"]
            kwargs["log_x"] = self.params["variable_options.log_values"]
        else:
            kwargs["log_y"] = self.params["variable_options.log_values"]

        # Set up the ordering
        category_orders = {
            cname: (
                data_long.groupby(cname)["value"].mean().sort_values().index
                if approach == "Mean"
                else (
                    data_long.groupby(cname)["value"].median().sort_values().index
                    if approach == "Median"
                    else data_long[cname].drop_duplicates().sort_values()
                )
            )
            for cname, approach in [
                (var_name, self.params["variable_options.sort_by"])
            ] + (
                [("category", self.params["category_options.sort_by"])]
                if category is not None
                else []
            )
        }

        fig = px.box(
            data_long,
            labels=dict(
                category=self.params["table.category.category.label"],
                value=self.params["display_options.val_label"],
                **{
                    var_name: self.params["display_options.var_label"]
                }
            ),
            points=(
                "outliers"
                if self.params["display_options.outliers.enabled"]
                else False
            ),
            category_orders=category_orders,
            **kwargs
        )
        fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
        fig.update_yaxes(matches=None)

        container.plotly_chart(fig)
