from typing import Optional
import pandas as pd
import plotly.express as px
from plotly.graph_objects import Figure
import streamlit as st
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
                    "enum": ["X-Axis", "Y-Axis", "Color", "Facet"],
                    "default": "Y-Axis",
                    "sidebar": True,
                },
                "log_values": {
                    "type": "boolean",
                    "label": "Log Scale",
                    "sidebar": True,
                },
                "sort_by": {
                    "type": "string",
                    "enum": ["Mean", "Median", "Name"],
                    "default": "Mean",
                    "label": "Sort By",
                    "sidebar": True,
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
                            "default": True,
                            "sidebar": True
                        }
                    },
                    "sidebar": True
                },
                "title": {
                    "type": "string",
                    "label": "Title",
                    "default": "",
                    "sidebar": True
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
                },
                "legend": {
                    "type": "string",
                    "label": "Legend",
                    "multiline": True
                }
            }
        }
    }
    legend_key = "display_options.legend"

    def build_figure(self) -> Figure:

        data: pd.DataFrame = self.params.get("table.data.dataframe")
        if data is None:
            st.write("Please select a data table")
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

        var_name = (
            "variable"
            if "variable" in data_long.columns
            else "var"
        )

        # If all of the columns start with the same prefix
        if data_long[var_name].apply(lambda x: x.split(":")[0]).nunique() == 1:
            data_long[var_name] = data_long[var_name].apply(
                lambda x: x.split(":", 1)[1] if ":" in x else x
            )

        # If the category was provided
        if category is not None:

            # Add it to the table
            data_long = data_long.assign(
                category=lambda df: df["index"].apply(category.get)
            )

        # Set up the kwargs
        kwargs = self._setup_kwargs(
            use_category=category is not None,
            var_name=var_name
        )
        if kwargs is None:
            return

        # Set up the ordering
        category_orders = {
            cname: (
                data_long.groupby(cname)["value"].mean().sort_values(ascending=False).index
                if approach == "Mean"
                else (
                    data_long.groupby(cname)["value"].median().sort_values(ascending=False).index
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
                category=self.params["table.category.columns.category.label"],
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
            color_discrete_sequence=px.colors.qualitative.D3,
            **kwargs
        )
        fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
        fig.update_yaxes(matches=None)

        return fig

    def _setup_kwargs(self, use_category: bool, var_name: str) -> Optional[dict]:

        # Set up the elements which are always the same
        kwargs = dict(
            facet_col_wrap=int(self.params["display_options.ncols"]),
            height=self.params["display_options.height"],
            title=self.params["display_options.title"]
        )

        # Make some more readable names for what we're working with
        category_axis = self.params["category_options.axis"] if use_category else None
        variable_axis = self.params["variable_options.axis"]
        log_values = self.params["variable_options.log_values"]

        # The x-axis will be the variable

        # If a category was provided
        if use_category:
            # If the category is on a facet
            if category_axis == "Facet":

                # Use the category as a facet
                kwargs['facet_col'] = "category"
                kwargs['boxmode'] = "overlay"

                # The variable cannot also be a facet
                if variable_axis == "Facet":
                    st.error("Cannot have both variable and category as facets.")
                    return
                
                # The variable cannot be a color, because then no axis is defined
                elif variable_axis == "Color":
                    st.error("Cannot have variable as color with category as facet - no axis is defined.")
                    return

                # If the  is on the x-axis, then the value is on the y-axis
                elif variable_axis == "X-Axis":
                    
                    kwargs["x"] = var_name
                    kwargs["y"] = "value"
                    kwargs["log_y"] = log_values

                # If the variable is on the y-axis, then the value is on the x-axis
                else:
                    if variable_axis != "Y-Axis":
                        st.error(f"Unexpected condition encountered (variable_axis == '{variable_axis}')")
                        return

                    kwargs["y"] = var_name
                    kwargs["x"] = "value"
                    kwargs["log_x"] = log_values

            # If the category is being shown on an axis
            elif category_axis == "Axis":

                # The variable cannot be on the x or y axis
                if variable_axis in ["X-Axis", "Y-Axis"]:
                    st.error("Cannot have variable as axis with category as axis.")
                    return
                
                # If the variable is a color
                elif variable_axis == "Color":
                    kwargs['x'] = "category"
                    kwargs['color'] = var_name
                    kwargs['y'] = 'value'
                    kwargs['log_y'] = log_values

                # If the variable is a facet and the category is on an axis
                else:
                    if variable_axis != "Facet":
                        st.error(f"Unexpected condition encountered (variable_axis == '{variable_axis}')")
                        return

                    kwargs['x'] = "category"
                    kwargs['y'] = 'value'
                    kwargs['log_y'] = log_values
                    kwargs['facet_col'] = var_name
                    kwargs['boxmode'] = "overlay"

            else:
                # Display the category as a color
                if category_axis != "Color":
                    st.error(f"Unexpected condition encountered (category_axis == '{category_axis}')")
                    return

                # The color is the category 
                kwargs['color'] = "category"

                # The variable cannot be a color
                if variable_axis == "Color":
                    st.error("Cannot have variable as color with category as color.")
                    return

                # The variable cannot be a facet
                elif variable_axis == "Facet":
                    st.error("Cannot have variable as facet with category as color.")
                    return

                # If the variable is on the x-axis, then the value is on the y-axis
                elif variable_axis == "X-Axis":
                    kwargs['x'] = var_name
                    kwargs['y'] = 'value'
                    kwargs['log_y'] = log_values

                # If the variable is on the y-axis, then the value is on the x-axis
                elif variable_axis == "Y-Axis":
                    kwargs['y'] = var_name
                    kwargs['x'] = 'value'
                    kwargs['log_x'] = log_values

                else:
                    st.error(f"Unexpected condition encountered (variable_axis == '{variable_axis}')")
                    return
                
        # If there is no category
        else:

            # If the variable is on the x-axis
            if variable_axis == "X-Axis":
                kwargs["x"] = var_name
                kwargs["y"] = "value"
                kwargs["log_y"] = log_values

            # If the variable is on the y-axis
            elif variable_axis == "Y-Axis":
                kwargs["y"] = var_name
                kwargs["x"] = "value"
                kwargs["log_x"] = log_values

            # The variable cannot be a color
            elif variable_axis == "Color":
                st.error("Cannot have variable as color with no category.")
                return
            
            # The variable cannot be a facet
            elif variable_axis == "Facet":
                st.error("Cannot have variable as facet with no category.")
                return

            else:
                st.error(f"Unexpected condition encountered (variable_axis == '{variable_axis}')")
                return
            
        return kwargs
