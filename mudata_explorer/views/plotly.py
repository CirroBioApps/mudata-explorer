import pandas as pd
from mudata_explorer.base.view import View
import plotly.express as px
from streamlit.delta_generator import DeltaGenerator


class Plotly(View):

    category = "Plotting"

    def fetch_dataframe(self, table_kw):
        """
        Helper method to fetch a dataframe from the params.
        This provides added functionality for parsing the color scales.
        """

        data: pd.DataFrame = self.params[f"{table_kw}.dataframe"]
        colorscale = {}

        # Check each column to see if it is a color scale
        for col_kw in [
            kw[:-len(".scale")] for kw in self.params.keys()
            if kw.endswith(".scale") and kw.startswith(table_kw)
        ]:
            cname = col_kw[len(table_kw) + 1:]

            if self.params[f"{col_kw}.enabled"]:
                if self.params[f"{col_kw}.is_categorical"]:
                    data = data.assign(**{
                        cname: data[cname].apply(str)
                    })
                    colorscale = dict(
                        color_discrete_sequence=getattr(
                            px.colors.qualitative,
                            self.params[f"{col_kw}.scale"]
                        )
                    )
                else:
                    colorscale = dict(
                        color_continuous_scale=self.params[f"{col_kw}.scale"]
                    )
            else:
                colorscale = {}

        return data, colorscale


class PlotlyScatter(Plotly):

    type = "plotly-scatter"
    name = "Scatterplot (Plotly)"
    help_text = "Display a two dimensional distribution of data using Plotly."
    schema = {
        "data": {
            "type": "dataframe",
            "columns": {
                "x": {"label": "x-axis"},
                "y": {"label": "y-axis"},
                "size": {"label": "size", "optional": True},
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

    def display(self, container: DeltaGenerator):

        data, colorscale = self.fetch_dataframe("data")

        fig = px.scatter(
            data.reset_index(),
            x="x",
            y="y",
            hover_name=(
                "index" if data.index.name is None else
                data.index.name
            ),
            log_x=self.params["scale_options.log_x"],
            log_y=self.params["scale_options.log_y"],
            color="color" if self.params["data.color.enabled"] else None,
            size="size" if self.params["data.size.enabled"] else None,
            labels=dict(
                x=self.params["data.x.label"],
                y=self.params["data.y.label"],
                color=self.params["data.color.label"],
                size=self.params["data.size.label"]
            ),
            **colorscale
        )

        container.plotly_chart(fig)


class PlotlyScatter3D(Plotly):

    type = "plotly-scatter-3d"
    name = "Scatterplot 3D (Plotly)"
    help_text = "Display a three dimensional distribution of data using Plotly."
    schema = {
        "data": {
            "type": "dataframe",
            "columns": {
                "x": {"label": "x-axis"},
                "y": {"label": "y-axis"},
                "z": {"label": "z-axis"},
                "size": {"label": "size", "optional": True},
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
                },
                "log_z": {
                    "type": "boolean",
                    "label": "Log Scale - Z Axis"
                }
            }
        }
    }

    def display(self, container: DeltaGenerator):

        data, colorscale = self.fetch_dataframe("data")

        fig = px.scatter_3d(
            data.reset_index(),
            x="x",
            y="y",
            z="z",
            hover_name=(
                "index" if data.index.name is None else
                data.index.name
            ),
            log_x=self.params["scale_options.log_x"],
            log_y=self.params["scale_options.log_y"],
            log_z=self.params["scale_options.log_z"],
            color="color" if self.params["data.color.enabled"] else None,
            size="size" if self.params["data.size.enabled"] else None,
            labels=dict(
                x=self.params["data.x.label"],
                y=self.params["data.y.label"],
                z=self.params["data.z.label"],
                color=self.params["data.color.label"],
                size=self.params["data.size.label"]
            ),
            **colorscale
        )

        container.plotly_chart(fig)


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
        }
    }

    def display(self, container: DeltaGenerator):

        data, colorscale = self.fetch_dataframe("data")

        if "color_continuous_scale" in colorscale:
            container.error("Color scale must be categorical for box plots.")

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
            **colorscale
        )

        container.plotly_chart(fig)


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
                }
            }
        }
    }

    def display(self, container: DeltaGenerator):

        data: pd.DataFrame = self.params["table.data.dataframe"]
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
            facet_col="variable",
            boxmode="overlay",
            facet_col_wrap=int(self.params["display_options.ncols"]),
            log_y=self.params["scale_options.log_y"],
            labels=dict(
                category=self.params["table.category.category.label"]
            )
        )
        fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
        fig.update_yaxes(matches=None)

        container.plotly_chart(fig)


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
            "help_text": "If a color column is used, should the bars be stacked or grouped?",
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
        }
    }

    def display(self, container: DeltaGenerator):

        data, colorscale = self.fetch_dataframe("data")

        if "color_continuous_scale" in colorscale:
            container.error("Color scale must be categorical for this display.")

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
            **colorscale
        )

        container.plotly_chart(fig)
        

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
                }
            }
        }
    }

    def display(self, container: DeltaGenerator):

        data: pd.DataFrame = self.params["table.data.dataframe"]
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
                        (c - c.min()) / (c.max() - c.min()) + 0.1
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

        fig = px.scatter(
            summary.reset_index(),
            x="column",
            y="category",
            size=size,
            color=color,
            hover_data=["Mean", "Median", "Positive", "Non-Null"],
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
