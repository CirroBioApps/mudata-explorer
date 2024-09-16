import os
import pandas as pd
from mudata_explorer.base.view import View
import plotly.express as px
from plotly.graph_objects import Figure
import streamlit as st


class Plotly(View):

    category = "Plotting"
    legend_key = "formatting.legend"

    def fetch_dataframe(self, table_kw):
        """
        Helper method to fetch a dataframe from the params.
        This provides added functionality for parsing the color scales.
        """

        if self.params.get(f"{table_kw}.dataframe") is None:
            return None, {}

        data: pd.DataFrame = self.params[f"{table_kw}.dataframe"]
        colorscale = {}

        # Check each column to see if it is a color scale
        for col_kw in [
            kw[:-len(".scale")] for kw in self.params.keys()
            if kw.endswith(".scale") and kw.startswith(table_kw + ".columns.")
        ]:
            cname = col_kw[len(table_kw) + len(".columns."):]

            if self.params.get(f"{col_kw}.enabled.value", True):
                if self.params[f"{col_kw}.is_categorical.value"]:
                    assert self.params[f"{col_kw}.scale.value"] is not None, \
                        f"Must specify color scale for {col_kw}"
                    data = data.assign(**{
                        cname: data[cname].apply(str)
                    })
                    colorscale = dict(
                        color_discrete_sequence=getattr(
                            px.colors.qualitative,
                            self.params[f"{col_kw}.scale.value"]
                        )
                    )
                else:
                    colorscale = dict(
                        color_continuous_scale=self.params[f"{col_kw}.scale.value"]
                    )
            else:
                colorscale = {}

        return data, colorscale
    
    def show_legend(self):
        """
        Helper method to show the legend.
        """
        legend = self.params.get(self.legend_key)
        if legend:
            st.markdown(legend)

    def build_figure(self) -> Figure:
        """
        By default, every plotly view will build a plotly figure.
        """

        raise NotImplementedError("Must implement build_figure method")

    def display(self):
        """
        By default, every plotly display will build a plotly chart and optionally
        display a legend.
        """

        fig = self.build_figure()
        if fig is not None:
            st.plotly_chart(fig)
            self.show_legend()

    def write_image(self, filename: str, **kwargs):
        """Write the image to a file."""
        fig = self.build_figure()
        assert fig is not None, "Figure is None"
        if fig is not None:
            fig.write_image(filename, **kwargs)
        assert os.path.exists(filename), f"File not found: {filename}"
