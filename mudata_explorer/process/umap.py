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
        mdata, modality, df, columns, use_zscore = inputs

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
                metric=metric
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

            # Add the UMAP coordinates to the obsm slot
            mdata.mod[modality].obsm[dest_key] = umap_df

            # Make a record of the process
            event = dict(
                process=self.type,
                params=params,
                timestamp=app.get_timestamp(),
                updated_keys=dest_key
            )

            # Update the MuData object
            app.set_mdata(mdata)

            # Add it to the history
            app.add_history(event)

            # Mark the source of the table which was added
            app.add_provenance(
                modality,
                "obsm",
                dest_key,
                event
            )

        # If the key already exists
        if dest_key in mdata.mod[modality].obsm.keys():
            prov = app.query_provenance(modality, "obsm", dest_key)
            if prov is not None:
                container.write(f"**Data currently in '{dest_key}'**")
                container.write(prov)
                container.write(
                    f"> 'Run' will overwrite existing data in '{dest_key}'."
                )


@st.cache_data
def run_umap(df: pd.DataFrame, **kwargs) -> pd.DataFrame:
    reducer = umap.UMAP(**kwargs)
    return pd.DataFrame(
        reducer.fit_transform(df),
        index=df.index,
        columns=["UMAP1", "UMAP2"]
    )
