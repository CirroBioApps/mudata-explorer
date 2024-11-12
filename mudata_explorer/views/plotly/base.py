from copy import copy
import json
import os
from typing import List
import pandas as pd
from mudata_explorer.base.view import View, ViewJSON
from mudata_explorer.base.form import _save_value
import plotly.express as px
from plotly.graph_objects import Figure
from scipy.stats import f_oneway, kruskal
import streamlit as st
from datastory.datastory import DataStory


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

    def _make_fig_list(self) -> List[Figure]:
        """
        Build figures and return as a list.
        """

        fig_list = self.build_figure()
        if fig_list is None:
            return
        
        if not isinstance(fig_list, list):
            fig_list = [fig_list]

        return fig_list

    def display(self):
        """
        By default, every plotly display will build a plotly chart and optionally
        display a legend.
        """

        # Get the list of figures
        fig_list = self._make_fig_list()

        # Display in a set of columns
        cols = st.columns(len(fig_list))

        for fig_ix, fig in enumerate(fig_list):
            # The selection will trigger a rerun only
            # if the form contains a MuDataAppPlotlySelection element
            on_select = (
                "rerun"
                if self.form.uses_plotly_selection(fig_ix)
                else "ignore"
            )
            # The input element can limit the selection modes available
            selection_mode = self.form.get_selection_mode(fig_ix)
            # Call the figure
            selection: dict = (
                cols[fig_ix]
                .plotly_chart(
                    copy(fig),
                    on_select=on_select,
                    selection_mode=selection_mode
                )
            )
            if self.form.uses_plotly_selection(fig_ix):
                if self.form.save_selection(selection.get("selection"), fig_ix):
                    self.save_changes()
                    st.rerun()

        self.show_legend()

    def to_json(self) -> ViewJSON:
        """Return a serialization of the displayed figure."""

        self.params = self.form.dump()

        fig_list = self._make_fig_list()

        if fig_list is None:
            return dict()

        return dict(
            type="plotly",
            figures=[
                json.loads(fig.to_json())
                for fig in fig_list
            ]
        )

    def write_image(self, filename: str, **kwargs):
        """Write the image to a file."""
        fig = self.build_figure()
        assert fig is not None, "Figure is None"
        if fig is not None:
            fig.write_image(filename, **kwargs)
        assert os.path.exists(filename), f"File not found: {filename}"

    def _add_stats_title(self, title: str, method: str, vals: list) -> str:
        """
        Add the results of a statistical test to the title of a figure.
        """
        if method == "ANOVA":
            try:
                res = f_oneway(*vals)
            except ValueError:
                return title
        elif method == "Kruskal-Wallis":
            try:
                res = kruskal(*vals)
            except ValueError:
                return title
        else:
            raise ValueError(f"Unknown method: {method}")

        # Format the p-value as a string
        pvalue = (
            f"{res.pvalue:.4f}"
            if res.pvalue > 0.0001
            else f"{res.pvalue:.2e}"
        )

        # Add it to the title
        formatted_result = f"{method} p-value: {pvalue}"
        if len(title) > 0:
            title = f"{title} ({formatted_result})"
        else:
            title = formatted_result

        return title

    def to_datastory(self, ds: DataStory):
        """
        Convert the view to a DataStory object.
        """
        
        self.params = self.form.dump()

        fig_list = self._make_fig_list()

        print(f"Number of figures for {self.name} #{self.ix +1}: {len(fig_list):,}")

        if len(fig_list) > 0:
            ds.add_section()
            ix = len(ds.sections) - 1

            for fig in fig_list:
                ds.add_plotly(fig, section_ix=ix)

        legend = self.params.get(self.legend_key)
        if legend:
            ds.add_markdown(legend, section_ix=ix)
