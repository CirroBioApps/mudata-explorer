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

    def run(self, container: DeltaGenerator):

        inputs = self.prompt_input_df(container)
        if inputs is None:
            return
        mdata, modality, axis, df, columns, use_zscore = inputs

        if container.checkbox(
            "Preview Clusters - Silhouette Scores",
            help="Evaluate a range of values for K using the silhouette score"
        ):
            # Try a number of different values for k
            self.show_silhouette_scores(df, container)

        # Prompt for the value of K to use
        k = container.number_input(
            "Number of Clusters (K)",
            min_value=2,
            value=5,
            step=1,
            help="Number of clusters to group samples into"
        )

        # Set the name of the obsm slot to use for the UMAP coordinates
        dest_key = container.text_input(
            "Destination Key",
            help="The name of the column which will be used for the results.",
            value="cluster"
        )

        # If the user clicks a button
        if container.button("Run K-Means Clustering"):

            # Run KMeans
            clusters = run_clustering(df, n_clusters=k)

            params = dict(
                dest_key=dest_key,
                modality=modality,
                axis=axis,
                columns=columns,
                use_zscore=use_zscore,
                k=k
            )

            # Save the results to the MuData object
            app.save_annot(
                mdata,
                modality,
                axis,
                dest_key,
                clusters,
                params,
                self.type
            )

        # Report to the user if data already exists in the destination key
        app.show_provenance(mdata, modality, axis, dest_key, container)

    def show_silhouette_scores(
        self,
        df: pd.DataFrame,
        container: DeltaGenerator
    ):
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
            n: run_clustering(df, n_clusters=n)
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
        container.plotly_chart(fig)


@st.cache_data
def run_clustering(df: pd.DataFrame, **kwargs):
    clusters = (
        KMeans(**kwargs)
        .fit_predict(df.values)
    )
    return pd.Series(
        clusters,
        index=df.index,
        name="cluster"
    ).apply(
        str
    )
