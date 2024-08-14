import numpy as np
import pandas as pd
from mudata_explorer.base.view import View
from mudata_explorer.app.mdata import get_mdata
import streamlit as st
from streamlit.delta_generator import DeltaGenerator
from typing import List, Optional
import plotly.express as px
from plotly import graph_objects as go
from plotly.subplots import make_subplots


class TSS_Tornado(View):

    type = "tss-tornado"
    name = "Transcription Start Site Browser"
    help_text = "Visualize chromatin accessibility surrounding transcription start sites."
    category = None
    schema = {}

    def runtime_options(self, container: DeltaGenerator):
        # Ask the user which display to show
        self._ask_view_type(container)

        # Ask the user whether to show the grouped TSS Tornado plot,
        # and if so, which value of K to show
        self._ask_selected_k(container)

        # Ask the user which clusters to display
        self._ask_selected_clusters(container)

        # Ask the user which samples to display
        self._ask_samples(container)

    def display(self, container: DeltaGenerator):

        if self.selected_k is not None:
            if len(self.selected_clusters) == 0:
                container.write("Please select >= 1 cluster to display")
                return

        if not self.samples or len(self.samples) == 0:
            container.write("Please select >= 1 sample to display")
            return

        if self.view_type == "UMAP":
            self.display_umap(container)
        elif self.view_type == "Fingerprint":
            self.display_fingerprint(container)
        elif self.view_type == "Trace":
            self.display_trace(container)
        elif self.view_type == "Trace + Fingerprint":
            self.display_trace_fingerprint(container)

    def display_umap(self, container: DeltaGenerator):
        mdata = get_mdata()
        # Get the UMAP coordinates
        umap = (
            mdata
            .mod["binned_coverage"]
            .obsm["umap"]
            .merge(
                mdata.obs,
                right_index=True,
                left_index=True
            )
        )

        # Filter to the selected samples
        umap = umap[umap["sample"].isin(self.samples)]

        # If a cluster was selected
        if self.selected_k is not None:
            # Annotate the cluster on the umap DataFrame
            umap = umap.assign(
                cluster=mdata.mod["binned_coverage"].obsm["kmeans"][str(self.selected_k)]
            )
            umap = umap[umap["cluster"].isin(self.selected_clusters)]

            # Plot the UMAP with the clusters
            container.plotly_chart(
                px.scatter(
                    umap,
                    x="UMAP1",
                    y="UMAP2",
                    color="cluster",
                    title=f"UMAP of Genes (K={self.selected_k})",
                    category_orders={"cluster": sorted(self.selected_clusters)}
                )
            )
        else:
            # Plot the UMAP without the clusters, coloring by sample
            container.plotly_chart(
                px.scatter(
                    umap,
                    x="UMAP1",
                    y="UMAP2",
                    color="sample",
                    title="UMAP of Genes",
                    category_orders={"sample": sorted(self.samples)}
                )
            )

    def display_trace(self, container: DeltaGenerator):
        """
        Make a line plot of the average coverage for the selected data slice.
        """

        fig = px.line(
            self.avg_coverage(),
            x="index",
            y="value",
            color="variable",
            labels={
                "index": "Distance from TSS (bp)",
                "value": "Average Coverage",
                "variable": "Sample",
                "cluster": "Cluster"
            },
            title=(
                "Average Coverage"
                if self.selected_k is None
                else f"Average Coverage (K={self.selected_k})"
            ),
            facet_row="cluster" if self.selected_k is not None else None,
            height=(
                600
                if self.selected_k is None
                else 400 * len(self.selected_clusters)
            )
        )
        fig.update_yaxes(matches=None)
        container.plotly_chart(fig)

    def display_trace_fingerprint(self, container: DeltaGenerator):
        """
        Make a line plot of the average coverage for the selected data slice,
        with a fingerprint plot below.
        """

        avg_coverage = self.avg_coverage()
        # Get the binned abundances
        mdata = get_mdata()
        bins = mdata.mod["binned_coverage"].to_df()

        if self.selected_k is None:
            # Show all of the data
            self.plot_tornado(
                avg_coverage,
                bins,
                mdata.obs,
                container,
                title="Average Coverage"
            )
        else:
            # Iterate over each of the selected clusters
            for cluster in self.selected_clusters:
                ix = mdata.mod["binned_coverage"].obsm["kmeans"][str(self.selected_k)] == cluster
                cluster_bins = bins.loc[ix]
                self.plot_tornado(
                    avg_coverage.query(f"cluster == '{cluster}'"),
                    cluster_bins,
                    mdata.obs.loc[ix],
                    container,
                    title=f"Average Coverage (K={self.selected_k}, {cluster})"
                )

    def plot_tornado(self, avg_coverage, bins, obs, container, title=""):
        fig = make_subplots(
            rows=2,
            shared_xaxes=True,
            shared_yaxes=False,
            row_heights=[0.3, 0.7]
        )
        for trace in px.line(
            avg_coverage,
            x="index",
            y="value",
            color="variable",
            labels={
                "index": "Distance from TSS (bp)",
                "value": "Average Coverage",
                "variable": "Sample",
                "cluster": "Cluster"
            }
        )["data"]:
            fig.add_trace(trace)

        fig.add_trace(
            go.Heatmap(
                z=bins.values,
                x=self.bin_position_labels,
                y=obs.apply(
                    lambda x: f"{x['name']} ({x['sample']})",
                    axis=1
                ),
                colorscale="blues",
                showscale=False
            ),
            row=2,
            col=1
        )
        fig.update_layout(
            title_text=title
        )

        container.plotly_chart(fig)

    def display_fingerprint(self, container: DeltaGenerator):
        # Get the binned abundances
        mdata = get_mdata()
        bins = mdata.mod["binned_coverage"].to_df()

        if self.selected_k is None:
            # Show all of the data
            self.plot_fingerprint(
                bins,
                mdata.obs,
                container
            )
        else:
            # Iterate over each of the selected clusters
            for cluster in self.selected_clusters:
                ix = mdata.mod["binned_coverage"].obsm["kmeans"][str(self.selected_k)] == cluster
                cluster_bins = bins.loc[ix]
                self.plot_fingerprint(
                    cluster_bins,
                    mdata.obs.loc[ix],
                    container,
                    title=f"K={self.selected_k}, {cluster}"
                )

    def plot_fingerprint(self, bins, obs, container, title=""):
        fig = make_subplots()

        fig.add_trace(
            go.Heatmap(
                z=bins.values,
                x=self.bin_position_labels,
                y=obs.apply(
                    lambda x: f"{x['name']} ({x['sample']})",
                    axis=1
                ),
                colorscale="blues",
                showscale=False
            ),
            row=1,
            col=1
        )
        fig.update_layout(title_text=title)

        container.plotly_chart(fig)

    @property
    def bin_position_labels(self):
        """Return the positions of the bins."""
        mdata = get_mdata()
        window_size = mdata.uns["avg_coverage"].shape[0]
        half_window = window_size // 2
        return np.arange(
            -half_window,
            half_window,
            window_size // mdata.mod["binned_coverage"].shape[1]
        )

    def avg_coverage(self):
        """Return the average coverage for the selected data slice."""

        # Get the mdata object
        mdata = get_mdata()

        # If we're showing for unclustered data
        if self.selected_k is None:
            # Get the per-sample means
            avg_coverage = (
                mdata
                .uns["avg_coverage"]
                .reindex(columns=self.samples)
                .fillna(0)
                .reset_index()
                .melt(id_vars=["index"])
            )

        else:
            avg_coverage = pd.concat(
                [
                    (
                        mdata
                        .uns[f"avg_coverage.{self.selected_k}.{selected_cluster}"]
                        .reindex(columns=self.samples)
                        .fillna(0)
                        .reset_index()
                        .melt(id_vars=["index"])
                        .assign(cluster=selected_cluster)
                    )                    
                    for selected_cluster in self.selected_clusters
                ]
            )

        return avg_coverage

    @property
    def selected_k(self) -> Optional[int]:
        val = self.uns.get("selected_k")
        if val is not None:
            return int(val)

    @selected_k.setter
    def selected_k(self, value: Optional[int]):
        if value != self.selected_k:
            self.update_view_uns("selected_k", value)
            self.selected_clusters = self.all_clusters

    def _ask_selected_k(self, params: DeltaGenerator):

        # The user can show the TSS Tornado plot for different values of K
        params.selectbox(
            "Grouping",
            ["All"] + [f"K={k}" for k in self.all_k],
            index=(
                0 if self.selected_k is None else self.selected_k - 1
            ),
            key=self.param_key("selected_k"),
            on_change=self._ask_selected_k_callback
        )

    def _ask_selected_k_callback(self):

        resp = st.session_state[self.param_key("selected_k")]

        if resp == "All":
            self.selected_k = None
        else:
            self.selected_k = int(resp[2:])

    @property
    def selected_clusters(self) -> Optional[List[str]]:
        # If no clustering method was selected, return None
        if self.selected_k is None:
            return

        return self.uns.get("selected_clusters", [])

    @selected_clusters.setter
    def selected_clusters(self, value: List[str]):
        self.update_view_uns("selected_clusters", value)

    def _ask_selected_clusters(self, params: DeltaGenerator):
        """
        Ask the user which of the clusters to display.
        """

        # If no clustering was selected, stop
        if self.selected_k is None:
            return

        # Make sure that the selection draws from the available
        # options
        self.selected_clusters = [
            c for c in self.selected_clusters if c in self.all_clusters
        ]

        params.multiselect(
            "Clusters",
            self.all_clusters,
            default=self.selected_clusters,
            key=self.param_key("selected_clusters"),
            on_change=self._ask_selected_clusters_callback
        )

    def _ask_selected_clusters_callback(self):
        self.selected_clusters = st.session_state[
            self.param_key("selected_clusters")
        ]

    @property
    def all_clusters(self) -> List[str]:
        if self.selected_k is None:
            return []
        # Get the clusters which are available
        options = (
            get_mdata()
            .mod["binned_coverage"]
            .obsm["kmeans"]
            [str(self.selected_k)]
            .drop_duplicates()
            .sort_values()
            .tolist()
        )
        options.sort()
        return options

    @property
    def view_type(self) -> str:
        return self.uns.get("view_type", "Trace")

    @view_type.setter
    def view_type(self, value: str):
        self.update_view_uns("view_type", value)

    def _ask_view_type(self, params: DeltaGenerator):
        """Ask the user what type of display to show."""
        options = [
            "Trace",
            "Trace + Fingerprint",
            "Fingerprint",
            "UMAP"
        ]

        # Get the display option which is saved
        if self.view_type not in options:
            self.view_type = options[0]

        params.selectbox(
            "View Type",
            options,
            index=options.index(self.view_type),
            key=self.param_key("view_type"),
            on_change=self._ask_view_type_callback
        )

    def _ask_view_type_callback(self):
        self.view_type = st.session_state[self.param_key("view_type")]

    @property
    def samples(self) -> List[str]:
        return [
            s
            for s in self.uns.get("samples", self.all_samples)
            if s in self.all_samples
        ]

    @samples.setter
    def samples(self, value: List[str]):
        self.update_view_uns("samples", value)

    def _ask_samples(self, params: DeltaGenerator):
        """Ask the user which samples to show."""

        params.multiselect(
            "Show Samples",
            self.all_samples,
            default=self.samples,
            key=self.param_key("samples"),
            on_change=self._ask_samples_callback
        )

    def _ask_samples_callback(self):
        self.samples = st.session_state[self.param_key("samples")]

    @property
    def all_samples(self) -> List[str]:
        # Get the MuData object
        mdata = get_mdata()
        if mdata is None:
            return []

        # Get all the samples which are available
        _all_samples = list(mdata.obs["sample"].unique())
        _all_samples.sort()
        return _all_samples

    @property
    def all_k(self) -> List[int]:
        # Get the MuData object
        mdata = get_mdata()
        if mdata is None:
            return []

        # Get all the values of K which are available
        all_k = [
            int(k)
            for k in (
                mdata.mod["binned_coverage"]
                .obsm["kmeans"]
                .columns
                .values
            )
        ]

        if len(all_k) == 0:
            return []
        all_k.sort()
        return all_k

    def _validate_mdata(self, mdata, container):
        if mdata is None:
            container.write("No data available.")
            return False

        for kw in ["chrom", "tss"]:
            if kw not in mdata.obs:
                container.write(f"Missing {kw} in metadata table.")
                return False

        # Make sure that the group means are available for all values of K
        for k in self.all_k:
            if f"kmeans_{k}_avg_zscore" not in mdata.uns:
                container.write(f"Missing group means for K={k}.")
                return False

        # Make sure that the global means are available
        if "avg_zscore" not in mdata.uns:
            container.write("Missing global means.")
            return False

        # Make sure that the selection of K is available
        selected_k = self.uns.get("selected_k")
        if selected_k is not None:
            if selected_k not in self.all_k:
                container.write(f"K={selected_k} not available.")
                return False

        return True
