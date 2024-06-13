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
            hover_name=data.index.name,
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
            hover_name=data.index.name,
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
    name = "Box Plot (Plotly)"
    help_text = "Display a series of data as a box graph using Plotly."
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
