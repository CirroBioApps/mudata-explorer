import json
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
        mdata = app.get_mdata()

        if mdata is None or mdata.shape[0] == 0:
            container.write("No MuData object available.")
            return

        # Select the modality to use
        modality = container.selectbox(
            "Select modality",
            list(mdata.mod.keys())
        )

        # Get the data for the selected modality
        df: pd.DataFrame = mdata.mod[modality].to_df()

        # If there is no data, return
        if df.shape[0] == 0 or df.shape[1] == 0:
            container.write(f"No data available for {modality}.")
            return

        # Let the user select the columns to use
        all_columns = list(df.columns.values)
        if container.checkbox("Use all columns", value=True):
            columns = all_columns
        else:
            columns = container.multiselect(
                "Select columns",
                all_columns,
                default=all_columns
            )

        if len(columns) == 0:
            container.write("No columns selected.")
            return

        df = df[columns].dropna()

        # Display the number of rows which contain values for all of the selected columns
        n_rows = df.shape[0]
        container.write(f"{n_rows:,} rows with data for all selected columns.")

        if n_rows == 0:
            return

        # Let the user optionally filter samples
        query = container.text_input(
            "Filter samples (optional)",
            help="Enter a query to filter samples (using metadata or data)."
        )
        if query is not None and len(query) > 0:
            df = (
                df
                .merge(mdata.obs, left_index=True, right_index=True)
                .query(query)
                .reindex(columns=columns)
                .dropna()
            )
            container.write(
                f"Filtered data: {df.shape[0]:,} rows x {df.shape[1]:,} columns."
            )

        if df.shape[0] == 0:
            return

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
        umap_key = container.text_input(
            "UMAP Key",
            value="X_umap"
        )
        # If the key already exists
        if umap_key in mdata.mod[modality].obsm.keys():
            provenance = app.query_provenance(modality, "obsm", umap_key)
            if provenance is not None:
                container.write(json.dumps(provenance, indent=4))
            else:
                container.write(
                    f"Warning: overwriting existing key '{umap_key}'."
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
                df[columns],
                **params
            )

            # Add the complete set of params
            params["umap_key"] = umap_key
            params["modality"] = modality
            params["columns"] = columns

            # Add the UMAP coordinates to the obsm slot
            mdata.mod[modality].obsm[umap_key] = umap_df

            # Make a record of the process
            event = dict(
                process=self.type,
                params=params,
                timestamp=str(pd.Timestamp.now()),
                updated_keys=umap_key
            )

            # Add it to the history
            app.add_history(event, mdata)

            # Mark the source of the table which was added
            app.add_provenance(
                umap_key,
                event,
                mdata
            )

            # Update the MuData object
            app.set_mdata(mdata)


@st.cache_data
def run_umap(df: pd.DataFrame, **kwargs) -> pd.DataFrame:
    reducer = umap.UMAP(**kwargs)
    return pd.DataFrame(
        reducer.fit_transform(df),
        index=df.index,
        columns=["UMAP1", "UMAP2"]
    )
