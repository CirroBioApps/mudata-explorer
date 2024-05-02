from sklearn.cluster import DBSCAN
import pandas as pd
import streamlit as st
from streamlit.delta_generator import DeltaGenerator
from mudata_explorer import app
from mudata_explorer.base.process import Process


class RunDBSCAN(Process):

    type = "dbscan"
    name = "DBSCAN"
    desc = "DBSCAN Clustering"
    categories = ["Clustering"]

    def run(self, container: DeltaGenerator):

        inputs = self.prompt_input_df(container)
        if inputs is None:
            return
        mdata, modality, df, columns, use_zscore = inputs

        # Prompt for the parameters to provide to DBSCAN
        eps = container.number_input(
            "Distance threshold - eps",
            min_value=0.,
            value=0.5,
            help="""
            The maximum distance between two samples for one to be considered
            as in the neighborhood of the other. This is not a maximum bound
            on the distances of points within a cluster. This is the most
            important DBSCAN parameter to choose appropriately for your data set
            and distance function.
            """
        )

        # Minimum number of samples
        min_samples = container.number_input(
            "Neighborhood Size - min_samples",
            min_value=2,
            step=1,
            value=5,
            help="""
            The number of samples (or total weight) in a neighborhood for a point to
            be considered as a core point. This includes the point itself. If
            `min_samples` is set to a higher value, DBSCAN will find denser clusters,
            whereas if it is set to a lower value, the found clusters will be more
            sparse.
            """
        )
        metric = container.selectbox(
            "DBSCAN: Metric",
            ["cosine", "euclidean", "manhattan", "correlation", "jaccard"],
            help="""
            The metric to use when calculating distance between instances in a
            feature array.
            """
        )

        # Set the name of the obsm slot to use for the UMAP coordinates
        dest_key = container.text_input(
            "Destination Key",
            help="The name of the column which will be used for the results.",
            value="cluster"
        )

        # If the user clicks a button
        if container.button("Run DBSCAN Clustering"):

            # Run KMeans
            clusters = run_clustering(
                df,
                eps=eps,
                min_samples=min_samples,
                metric=metric
            )

            # Add the complete set of params
            params = dict(
                dest_key=dest_key,
                modality=modality,
                columns=columns,
                use_zscore=use_zscore,
                eps=eps,
                min_samples=min_samples,
                metric=metric
            )

            # Add the PCA coordinates to the obsm slot
            mdata.mod[modality].obs[dest_key] = clusters

            # Report the number of clusters
            container.write("#### Number of observations per cluster")
            container.write(pd.Series(clusters).value_counts())

            # Make a record of the process
            event = dict(
                process=self.type,
                params=params,
                timestamp=app.get_timestamp(),
                updated_keys=dest_key
            )

            # Save the results

            # Update the MuData object
            app.set_mdata(mdata)

            # Add it to the history
            app.add_history(event)

            # Mark the source of the table which was added
            app.add_provenance(
                modality,
                "obs",
                dest_key,
                event
            )

        # If the key already exists
        if dest_key in mdata.mod[modality].obs.keys():
            prov = app.query_provenance(modality, "obs", dest_key)
            if prov is not None:
                container.write(f"**Data currently in '{dest_key}'**")
                container.write(prov)
                container.write(
                    f"> 'Run' will overwrite existing data in '{dest_key}'."
                )


@st.cache_data
def run_clustering(df: pd.DataFrame, **kwargs):
    clusters = (
        DBSCAN(**kwargs)
        .fit_predict(df.values)
    )
    return list(map(str, clusters))
