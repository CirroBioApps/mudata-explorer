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

    def run(self, container: DeltaGenerator):

        inputs = self.prompt_input_df(container)
        if inputs is None:
            return
        mdata, modality, axis, df, columns, use_zscore = inputs
        slot = axis + "m"

        n_neighbors = container.number_input(
            "UMAP: Number of neighbors",
            value=15
        )
        min_dist = container.number_input(
            "UMAP: Minimum distance",
            value=0.1
        )
        metric = container.selectbox(
            "UMAP: Metric",
            ["cosine", "euclidean", "manhattan", "correlation", "jaccard"]
        )
        n_components = container.number_input(
            "UMAP: Number of Components",
            value=2,
            min_value=1,
            step=1,
            help="The number of dimensions to reduce the data to.",
        )

        # Set the name of the obsm slot to use for the UMAP coordinates
        dest_key = container.text_input(
            "Destination Key",
            help="The name of the table which will be used for the results.",
            value="X_umap"
        )

        # If the user clicks a button
        if container.button("Run UMAP"):

            params = dict(
                n_neighbors=n_neighbors,
                min_dist=min_dist,
                metric=metric,
                n_components=n_components
            )

            # Run UMAP
            umap_df = run_umap(
                df,
                **params
            )

            # Add the complete set of params
            params["dest_key"] = dest_key
            params["modality"] = modality
            params["columns"] = columns
            params["use_zscore"] = use_zscore

            # Save the results to the MuData object
            app.save_annot(
                mdata,
                modality,
                slot,
                dest_key,
                umap_df,
                params,
                self.type
            )

        # Report to the user if data already exists in the destination key
        app.show_provenance(mdata, modality, slot, dest_key, container)


@st.cache_data
def run_umap(df: pd.DataFrame, n_components=2, **kwargs) -> pd.DataFrame:
    reducer = umap.UMAP(n_components=n_components, **kwargs)
    return pd.DataFrame(
        reducer.fit_transform(df),
        index=df.index,
        columns=[f"UMAP-{i+1}" for i in range(n_components)]
    )
