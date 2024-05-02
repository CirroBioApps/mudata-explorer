from typing import List, Tuple, Union
import pandas as pd
import muon as mu
from mudata_explorer import app
from scipy.stats import zscore
from streamlit.delta_generator import DeltaGenerator


class Process:

    type: str
    name: str
    desc: str
    categories: List[str]

    def run(self, container: DeltaGenerator):
        pass

    def prompt_input_df(
        self,
        container: DeltaGenerator
    ) -> Union[
        None, 
        Tuple[
            mu.MuData,
            str,
            pd.DataFrame,
            List[str]
        ]
    ]:
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

        use_zscore = container.checkbox(
            "Normalize: Z-score",
            value=False,
            help="Weight each variable equally by computing the z-score"
        )

        df = df[columns]

        if use_zscore:
            # Drop any columns with no variance
            dropping_cnames = [
                cname for cname in columns
                if df[cname].std() == 0
            ]
            if len(dropping_cnames) > 0:
                container.write("Dropping columns with zero variance:")
                container.write("- " + "\n- ".join(dropping_cnames))
                columns = list(set(columns) - set(dropping_cnames))

            df = df.apply(zscore)

        return mdata, modality, df, columns, use_zscore
