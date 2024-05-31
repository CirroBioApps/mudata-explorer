from typing import Union
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import pandas as pd
import plotly.express as px
from plotly import io
import streamlit as st
from streamlit.delta_generator import DeltaGenerator
from mudata_explorer import app
from mudata_explorer.base.process import Process


class RunKmeans(Process):

    type = "kmeans"
    name = "K-Means"
    desc = "K-Means Clustering"
    category = "Clustering"
    output_type = pd.Series
    schema = {
        "data": {
            "type": "dataframe",
            "select_columns": True,
            "query": "",
        },
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
            "label": "Evaluate Values of K: Lower Bound",
            "help": """
            Generate a plot showing the silhouette score for a range of K
            """
        },
        "max_k": {
            "type": "integer",
            "min_value": 2,
            "default": 10,
            "label": "Evaluate Values of K: Upper Bound",
            "help": """
            Generate a plot showing the silhouette score for a range of K
            """
        }
    }

    def execute(self) -> Union[pd.Series, pd.DataFrame]:

        msg = "Selected value of K must fall between the lower and upper bound"
        assert self.params["k"] >= self.params["min_k"], msg
        assert self.params["k"] <= self.params["max_k"], msg

        df: pd.DataFrame = self.params["data.dataframe"].dropna()
        msg = "Null values in all rows - remove invalid columns"
        assert df.shape[0] > 0, msg

        # Cluster the data
        clusters = {
            n: run_clustering(
                df,
                n_clusters=n
            )
            for n in range(
                self.params["min_k"],
                self.params["max_k"]
            )
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
        fig.add_vline(x=self.params["k"], line_dash="dash", line_color="grey")

        # Save the figure
        self.figures = [io.to_json(fig, validate=False)]

        # Return the clusters for the selected value of K
        return clusters[self.params["k"]]


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
