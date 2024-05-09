import pandas as pd
from mudata_explorer.base.view import View
import plotly.express as px
from streamlit.delta_generator import DeltaGenerator


class Plotly(View):

    categories = ["Plotting"]


class PlotlyScatter(Plotly):

    type = "plotly-scatter"
    name = "Scatterplot (Plotly)"
    desc = "Display a two dimensional distribution of data using Plotly."
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
                    "continuous_scale": True
                },
            },
            "query": "",
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

        # Parse information from the user,
        # as defined by the form schema
        self.get_data(container)

        data: pd.DataFrame = self.params["data.dataframe"]

        fig = px.scatter(
            data,
            x=self.params["data.x"],
            y=self.params["data.y"],
            log_x=self.params["scale_options.log_x"],
            log_y=self.params["scale_options.log_y"],
            color=self.params["data.color"] if self.params["data.color.enabled"] else None,
            color_continuous_scale=self.params["data.color.continuous_scale"],
            size=self.params["data.size"] if self.params["data.size.enabled"] else None
        )

        container.plotly_chart(fig)


class PlotlyScatter3D(Plotly):

    type = "plotly-scatter-3d"
    name = "Scatterplot 3D (Plotly)"
    desc = "Display a three dimensional distribution of data using Plotly."
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
                    "continuous_scale": True
                },
            },
            "query": "",
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

        # Parse information from the user,
        # as defined by the form schema
        self.get_data(container)

        data: pd.DataFrame = self.params["data.dataframe"]

        fig = px.scatter_3d(
            data,
            x=self.params["data.x"],
            y=self.params["data.y"],
            z=self.params["data.z"],
            log_x=self.params["scale_options.log_x"],
            log_y=self.params["scale_options.log_y"],
            log_z=self.params["scale_options.log_z"],
            color=self.params["data.color"] if self.params["data.color.enabled"] else None,
            color_continuous_scale=self.params["data.color.continuous_scale"],
            size=self.params["data.size"] if self.params["data.size.enabled"] else None
        )

        container.plotly_chart(fig)


class PlotlyLine(Plotly):

    type = "plotly-line"
    name = "Line Plot (Plotly)"
    desc = "Display a series of data as a line graph using Plotly."
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
                    "discrete_sequence": True
                },
            },
            "query": "",
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
        # Parse information from the user,
        # as defined by the form schema
        self.get_data(container)

        data: pd.DataFrame = self.params["data.dataframe"]
        data = data.sort_values(self.params["data.sort_by"])

        fig = px.line(
            data,
            x=self.params["data.x"],
            y=self.params["data.y"],
            color=self.params["data.color"],
            color_discrete_sequence=self.params["data.color.discrete_sequence"],
            log_x=self.params["scale_options.log_x"],
            log_y=self.params["scale_options.log_y"],
        )

        container.plotly_chart(fig)


class PlotlyBox(Plotly):

    type = "plotly-box"
    name = "Box Plot (Plotly)"
    desc = "Display a series of data as a box graph using Plotly."
    schema = {
        "data": {
            "type": "dataframe",
            "columns": {
                "x": {"label": "x-axis"},
                "y": {"label": "y-axis"},
                "color": {
                    "label": "Color",
                    "optional": True,
                    "discrete_sequence": True
                },
            },
            "query": "",
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

        # Parse information from the user,
        # as defined by the form schema
        self.get_data(container)

        data: pd.DataFrame = self.params["data.dataframe"]

        fig = px.box(
            data,
            x=self.params["x"],
            y=self.params["y"],
            color=self.params["color"],
            log_y=self.params["log_y"]
        )

        container.plotly_chart(fig)
