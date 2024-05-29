from typing import Union
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import pandas as pd
import plotly.express as px
import streamlit as st
from streamlit.delta_generator import DeltaGenerator
from mudata_explorer import app
from mudata_explorer.base.process import Process


class RunKmeans(Process):

    type = "kmeans"
    name = "K-Means"
    desc = "K-Means Clustering"
    categories = ["Clustering"]
    output_type = pd.Series
    schema = {
        "data": {
            "type": "dataframe",
            "select_columns": True,
            "query": "",
        },
        "preview_clusters": {
            "type": "boolean",
            "label": "Preview Clusters - Silhouette Scores",
            "help": """
            Evaluate a range of values for K using the silhouette score
            """
        },
        "k": {
            "type": "integer",
            "min_value": 2,
            "value": 5,
            "label": "Number of Clusters (K)",
            "help": """
            Number of clusters to group samples into
            """
        }
    }

    def execute(self) -> Union[pd.Series, pd.DataFrame]:

        # Run KMeans
        return run_clustering(
            self.params["data.dataframe"],
            n_clusters=self.params["k"]
        )

    def display(self, container: DeltaGenerator):
        if not self.params["preview_clusters"]:
            return

        min_k = container.number_input(
            "Min: K",
            min_value=2,
            value=2,
            step=1,
            help="Smallest value of K used for evaluation"
        )
        max_k = container.number_input(
            "Min: K",
            min_value=2,
            value=10,
            step=1,
            help="Largest value of K used for evaluation"
        )

        # Cluster the data
        clusters = {
            n: run_clustering(
                self.params["data.dataframe"],
                n_clusters=n
            )
            for n in range(min_k, max_k)
            if n < self.params["data.dataframe"].shape[0]
        }

        # Compute silhouette scores for the clusters
        silhouette_scores = {
            n: silhouette_score(self.params["data.dataframe"], clust)
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
        container.plotly_chart(fig)


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
