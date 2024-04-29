import pandas as pd
import umap
import streamlit as st
from streamlit.delta_generator import DeltaGenerator
from mudata_explorer.base.process import Process
from muon import MuData


class UMAP(Process):

    type = "umap"
    name = "UMAP"
    desc = "Uniform Manifold Approximation and Projection (UMAP)"
    categories = ["Dimensionality Reduction"]

    def run(self, container: DeltaGenerator, mdata: MuData):
        mdata = self.get_mdata()

        if mdata is None or mdata.shape[0] == 0:
            container.write("No MuData object available.")
            return

        return

        n_neighbors = container.number_input(
            "UMAP: Number of neighbors",
            value=15,
            key=self.param_key("n_neighbors")
        ),
        min_dist = container.number_input(
            "UMAP: Minimum distance",
            value=0.1,
            key=self.param_key("min_dist")
        ),
        metric = container.selectbox(
            "UMAP: Metric",
            ["cosine", "euclidean", "manhattan", "correlation", "jaccard"],
            key=self.param_key("metric")
        )



@st.cache_data
def run_umap(df: pd.DataFrame, **kwargs) -> pd.DataFrame:
    reducer = umap.UMAP(**kwargs)
    return pd.DataFrame(
        reducer.fit_transform(df),
        index=df.index,
        columns=["UMAP1", "UMAP2"]
    )
