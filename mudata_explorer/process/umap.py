import pandas as pd
import umap
import streamlit as st
from streamlit.delta_generator import DeltaGenerator
from mudata_explorer import app
from mudata_explorer.base.process import Process


class UMAP(Process):

    type = "umap"
    name = "UMAP"
    desc = "Uniform Manifold Approximation and Projection (UMAP)"
    categories = ["Dimensionality Reduction"]
    output_type = pd.DataFrame
    schema = {
        "data": {
            "type": "dataframe",
            "select_columns": True,
            "query": "",
        },
        "n_neighbors": {
            "type": "integer",
            "min_value": 2,
            "default": 15,
            "label": "Number of Neighbors",
            "help": "Number of neighbors to use in the UMAP algorithm."
        },
        "min_dist": {
            "type": "float",
            "min_value": 0.,
            "default": 0.1,
            "label": "Minimum Distance",
            "help": "Minimum distance between points in the UMAP algorithm."
        },
        "metric": {
            "type": "string",
            "label": "UMAP: Metric",
            "enum": ["cosine", "euclidean", "manhattan", "correlation", "jaccard"],
            "default": "correlation",
            "help": "The metric to use for the UMAP algorithm."
        },
        "n_components": {
            "type": "integer",
            "min_value": 1,
            "default": 2,
            "label": "Number of Components",
            "help": "Number of dimensions to reduce the data to."
        }
    }

    def execute(self) -> pd.DataFrame:

        df: pd.DataFrame = self.params["data.dataframe"]

        # Run UMAP
        return run_umap(
            df,
            **{
                kw: self.params[kw]
                for kw in [
                    "n_neighbors",
                    "min_dist",
                    "metric",
                    "n_components"
                ]
            }
        )


@st.cache_data
def run_umap(df: pd.DataFrame, n_components=2, **kwargs) -> pd.DataFrame:
    reducer = umap.UMAP(n_components=n_components, **kwargs)
    return pd.DataFrame(
        reducer.fit_transform(df),
        index=df.index,
        columns=[f"UMAP {i+1}" for i in range(n_components)]
    )
