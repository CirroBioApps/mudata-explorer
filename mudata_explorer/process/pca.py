from sklearn.decomposition import PCA
from scipy.stats import zscore
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
        mdata, modality, df, columns, use_zscore = inputs

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

            # Add the PCA coordinates to the obsm slot
            mdata.mod[modality].obsm[dest_key] = pca_df

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
def run_pca(df: pd.DataFrame) -> pd.DataFrame:
    pca = PCA()
    return pd.DataFrame(
        pca.fit_transform(df),
        index=df.index,
        columns=[f"PC{i + 1}" for i in range(df.shape[1])]
    )
