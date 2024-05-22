from sklearn.decomposition import PCA
import pandas as pd
import streamlit as st
from streamlit.delta_generator import DeltaGenerator
from mudata_explorer import app
from mudata_explorer.base.process import Process


class RunPCA(Process):

    type = "pca"
    name = "PCA"
    desc = "Principle Coordinates Analysis (PCA)"
    categories = ["Dimensionality Reduction"]

    def run(self, container: DeltaGenerator):

        inputs = self.prompt_input_df(container)
        if inputs is None:
            return
        mdata, modality, axis, df, columns, use_zscore = inputs

        slot = axis + "m"

        # Set the name of the obsm slot to use for the UMAP coordinates
        dest_key = container.text_input(
            "Destination Key",
            help="The name of the table which will be used for the results.",
            value="X_pca"
        )

        # If the user clicks a button
        if container.button("Run PCA"):

            # Run PCA
            pca_df = run_pca(df)

            # Add the complete set of params
            params = dict(
                dest_key=dest_key,
                modality=modality,
                columns=columns,
                use_zscore=use_zscore
            )

            # Save the results to the MuData object
            app.save_annot(
                mdata,
                modality,
                slot,
                dest_key,
                pca_df,
                params,
                self.type
            )

        # Report to the user if data already exists in the destination key
        app.show_provenance(mdata, modality, slot, dest_key, container)


@st.cache_data
def run_pca(df: pd.DataFrame) -> pd.DataFrame:
    pca = PCA()
    ndim = min(df.shape[0], df.shape[1])
    return pd.DataFrame(
        pca.fit_transform(df),
        index=df.index,
        columns=[f"PC{i + 1}" for i in range(ndim)]
    )
