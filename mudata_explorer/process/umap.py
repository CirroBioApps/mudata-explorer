import pandas as pd
import umap
import streamlit as st
from mudata_explorer.base.process import Process


class UMAP(Process):

    type = "umap"
    name = "Uniform Manifold Approximation and Projection (UMAP)"
    help_text = """
    UMAP is a fairly flexible non-linear dimension reduction algorithm.
    It seeks to learn the manifold structure of your data and find a low
    dimensional embedding that preserves the essential topological structure
    of that manifold.

    Compared to something like PCA, UMAP is able to capture more complex
    relationships between groups of data points using a smaller number of
    dimensions.
    Instead of having to look through PC1, PC2, PC3, etc. to understand the
    structure of your data, UMAP can often capture the same information in
    just UMAP1 and UMAP2.

    - [UMAP Documentation](https://umap-learn.readthedocs.io/en/latest/)
    - [Understanding UMAP](https://pair-code.github.io/understanding-umap/)
    - [UMAP in Single-Cell Analysis](https://alleninstitute.org/resource/what-is-a-umap/)
    """ # noqa
    category = "Dimensionality Reduction"
    schema = {
        "table": {
            "type": "object",
            "properties": {
                "data": {
                    "label": "Data Table",
                    "type": "dataframe",
                    "select_columns": True,
                    "query": "",
                }
            }
        },
        "umap_params": {
            "type": "object",
            "label": "UMAP Parameters",
            "properties": {
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
                    "help": "Minimum distance between points."
                },
                "metric": {
                    "type": "string",
                    "label": "UMAP: Metric",
                    "enum": [
                        "braycurtis",
                        "canberra",
                        "chebyshev",
                        "correlation",
                        "cosine",
                        "euclidean",
                        "haversine",
                        "mahalanobis",
                        "manhattan",
                        "minkowski",
                        "seuclidean",
                        "wminkowski"
                    ],
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
        },
        "outputs": {
            "type": "object",
            "label": "Outputs",
            "properties": {
                "dest_key": {
                    "type": "string",
                    "default": "X_umap",
                    "label": "Label to use for results",
                    "help": """
                    Key to use when saving the output
                    """
                }
            }
        }
    }
    outputs = {
        "res": {
            "type": pd.DataFrame,
            "label": "UMAP Coordinates",
            "desc": "Low dimensional embedding returned by UMAP",
            "modality": "table.data.tables",
            "axis": "table.data.axis",
            "attr": "outputs.dest_key"
        }
    }

    def execute(self):

        df: pd.DataFrame = self.params["table.data.dataframe"]
        assert df is not None, self.params

        # Run UMAP
        res = run_umap(
            df,
            **{
                kw: self.params[f"umap_params.{kw}"]
                for kw in [
                    "n_neighbors",
                    "min_dist",
                    "metric",
                    "n_components"
                ]
            }
        )
        # Save the results and the figure
        self.save_results("res", res)


@st.cache_data
def run_umap(df: pd.DataFrame, n_components=2, **kwargs) -> pd.DataFrame:
    reducer = umap.UMAP(n_components=n_components, **kwargs)
    return pd.DataFrame(
        reducer.fit_transform(df),
        index=df.index,
        columns=[f"UMAP {i+1}" for i in range(n_components)]
    )
