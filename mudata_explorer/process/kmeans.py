from typing import Union
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import pandas as pd
import plotly.express as px
from plotly import io
import streamlit as st
from mudata_explorer.base.process import Process


class RunKmeans(Process):

    type = "kmeans"
    name = "K-Means Clustering"
    help_text = """
    K-Means clustering is a method of vector quantization that aims to partition
    n observations into k clusters in which each observation belongs to the
    cluster with the nearest mean.

    The user selects the number of clusters (k) to group the samples into.

    Silhouette scores are computed over a range of values of k to evaluate the
    clustering performance. The silhouette score is a measure of how similar
    an object is to its own cluster compared to other clusters.

    - [Wikipedia: K-Means Clustering](https://en.wikipedia.org/wiki/K-means_clustering)
    - [Visual Explanation of K-Means](https://www.naftaliharris.com/blog/visualizing-k-means-clustering/)
    """ # noqa
    category = "Clustering"
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
        "clustering": {
            "type": "object",
            "label": "Clustering Parameters",
            "properties": {
                "k": {
                    "type": "integer",
                    "min_value": 2,
                    "default": 5,
                    "label": "Number of Clusters (K)",
                    "help": """
                    Number of clusters to group samples into
                    """
                },
                "min_k": {
                    "type": "integer",
                    "min_value": 2,
                    "default": 2,
                    "label": "Evaluate Values of K: Lower Bound"
                },
                "max_k": {
                    "type": "integer",
                    "min_value": 2,
                    "default": 10,
                    "label": "Evaluate Values of K: Upper Bound"
                }
            }
        },
        "outputs": {
            "type": "object",
            "label": "Outputs",
            "properties": {
                "dest_key": {
                    "type": "string",
                    "default": "kmeans",
                    "label": "Label to use for results",
                    "help": """
                    Key to use when saving the output to the container
                    """
                }
            }
        }
    }
    outputs = {
        "res": {
            "type": pd.Series,
            "label": "Clusters",
            "desc": "Cluster assignments for each element",
            "modality": "table.data.tables",
            "axis": "table.data.axis",
            "attr": "outputs.dest_key"
        }
    }

    def execute(self) -> Union[pd.Series, pd.DataFrame]:

        k = self.params["clustering.k"]
        min_k = self.params["clustering.min_k"]
        max_k = self.params["clustering.max_k"]

        msg = "Selected value of K must fall between the lower and upper bound"
        assert k >= min_k, msg
        assert k <= max_k, msg

        df: pd.DataFrame = self.params["table.data.dataframe"].dropna()
        msg = "Null values in all rows - remove invalid columns"
        assert df.shape[0] > 0, msg

        # Cluster the data
        clusters = {
            n: run_clustering(
                df,
                n_clusters=n
            )
            for n in range(min_k, max_k)
            if n < df.shape[0]
        }

        # Compute silhouette scores for the clusters
        silhouette_scores = {
            n: silhouette_score(df, clust)
            for n, clust in clusters.items()
        }

        # Plot the silhouette scores
        fig = px.line(
            x=list(silhouette_scores.keys()),
            y=list(silhouette_scores.values()),
            labels={
                "x": "Number of clusters",
                "y": "Silhouette score (1=best)"
            },
            title="Evaluate Clustering Performance"
        )
        fig.add_vline(x=k, line_dash="dash", line_color="grey")

        # Save the results and the figure
        self.save_results(
            "res",
            clusters[k],
            figures=[io.to_json(fig, validate=False)]
        )


@st.cache_data
def run_clustering(df: pd.DataFrame, **kwargs):
    clusters = (
        KMeans(**kwargs)
        .fit_predict(df.values)
    )
    return pd.Series(
        list(map(str, clusters)),
        index=df.index
    )
