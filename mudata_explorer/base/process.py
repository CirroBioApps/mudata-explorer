from typing import List, Tuple, Union
import pandas as pd
import muon as mu
from mudata_explorer import app
from mudata_explorer.base.base import MuDataAppHelpers
from scipy.stats import zscore
from streamlit.delta_generator import DeltaGenerator


class Process(MuDataAppHelpers):

    type: str
    name: str
    desc: str
    categories: List[str]
    schema: dict
    ix = -1

    def __init__(
        self,
        params={}
    ):
        self.params = {
            kw: params.get(kw, val)
            for kw, val in self.get_schema_defaults(self.schema)
        }

    def run(self, container: DeltaGenerator):

        # Get the parameters from the user
        self.get_data(container)

        container.write(self.params)

        # # Now run the method, catching any errors
        # try:
        #     self.display(self.view_container)
        # except Exception as e:
        #     # Log the full traceback of the exception
        #     self.view_container.exception(e)

    def param_key(self, kw):
        return f"process-{kw}"

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

        # Select the axis to use
        axis = container.selectbox(
            "Select axis",
            ["Observations", "Variables"],
            help="Select the axis to use for the analysis."
        )
        axis = axis.lower()[:3]
        if axis == "var":
            df = df.T

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

        # Display the number of rows which contain values
        # for all of the selected columns
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
                f"Filtered: {df.shape[0]:,} rows x {df.shape[1]:,} columns."
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

        return mdata, modality, axis, df, columns, use_zscore

    def update_view_param(self, kw, value):
        # Get the MuData object
        mdata = app.get_mdata()

        # Modify the value of this param for this view
        mdata.uns["mudata-explorer-process"]["params"][kw] = value

        # Save the MuData object
        app.set_mdata(mdata)

        # Also update the params object
        self.params[kw] = value
