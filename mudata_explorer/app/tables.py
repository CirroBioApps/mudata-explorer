import streamlit as st
from streamlit.delta_generator import DeltaGenerator
import pandas as pd
from mudata import MuData
from mudata_explorer.app.mdata import set_mdata
from anndata import AnnData


def _sanitize_types(
    df: pd.DataFrame,
    container: DeltaGenerator,
    keep_str=False
):
    # Convert every column to the type of the first non-null value
    to_drop = []
    for cname in df.columns:
        if df[cname].dtype == object:
            if isinstance(df[cname].dropna().values[0], str):
                if keep_str:
                    df[cname] = df[cname].fillna("").apply(str)
                    df.dtypes[cname] = str
                else:
                    to_drop.append(cname)
            else:
                try:
                    df[cname] = df[cname].apply(
                        type(df[cname].dropna().values[0])
                    )
                except Exception:
                    to_drop.append(cname)
        elif df[cname].dtype == str:
            to_drop.append(cname)
    if len(to_drop) > 0:
        container.write(f"Dropping non-numeric columns: {to_drop}")
        df = df.drop(columns=to_drop)

    return df


def get_table(
    container: DeltaGenerator,
    label: str,
    help=None,
    keep_str=False,
    key="_get_table"
):

    # Get input from the user
    file = container.file_uploader(
        label,
        help=help,
        key=key
    )
    df = _read_table(file, container)
    if df is not None:
        df = _sanitize_types(df, container, keep_str=keep_str)
    return df


def _read_table(file, container: DeltaGenerator):

    if file is None:
        return

    if file.name.endswith("xlsx"):
        try:
            df = pd.read_excel(file)
        except pd.errors.EmptyDataError:
            container.write("Could not read data from file")
            return
    elif file.name.endswith("csv"):
        try:
            df = pd.read_csv(file)
        except pd.errors.EmptyDataError:
            container.write("Could not read data from file")
            return
    elif file.name.endswith("tsv"):
        try:
            df = pd.read_csv(file, sep="\t")
        except pd.errors.EmptyDataError:
            container.write("Could not read data from file")
            return
    else:
        container.error("File must be a CSV or TSV.")
        return

    # The first column must only have unique values
    if not df.iloc[:, 0].is_unique:
        container.error("The first column must have unique values.")
        return

    # Make sure that the first column is interpreted as a string
    index_col = df.columns.values[0]
    df = df.assign(**{index_col: df[index_col].apply(str)})

    # Set the first column as the index
    df = df.set_index(index_col)

    # Return the data
    return df


def build_mdata(id="main"):

    # Let the user upload a metadata table
    obs_df = get_table(
        st.container(),
        "Upload Observation Metadata (CSV, TSV, or XLSX)",
        keep_str=True,
        key="obs"
    )
    if obs_df is not None:
        st.write(f"Read in {obs_df.shape[0]:,} rows x {obs_df.shape[1]:,} columns.")

    # Let the user specify the number of modalities to add
    n_mods = st.number_input(
        "Number of Measurement Modalities",
        min_value=1,
        value=1
    )

    mod_data = []

    for n in range(n_mods):
        mod_name = st.text_input(
            f"Measurement Name {n + 1}",
            help="Enter the name of the measurement.",
            value=f"Measurement {n + 1}"
        )
        if "." in mod_name:
            st.error("The measurement name cannot contain a period.")
            return
        mod_df = get_table(
            st.container(),
            f"Upload New Table {n + 1}",
            key=f"mod_{n}"
        )
        if mod_df is not None:
            mod_data.append((mod_name, mod_df))
            st.write(
                f"Read in {mod_df.shape[0]:,} rows x {mod_df.shape[1]:,} columns."
            )

    if obs_df is None:
        st.write("Note: No observation metadata was uploaded.")
        samples = set()
    else:
        samples = set(obs_df.index)

    # Get the total unique set of samples across all modalities
    for mod_name, mod_df in mod_data:
        samples.update(mod_df.index)
    for mod_name, mod_df in mod_data:
        st.write(f"{mod_name}: {mod_df.shape[0]:,} / {len(samples):,} total samples present")

    if obs_df is not None:
        obs_df = obs_df.reindex(index=samples)

    if len(mod_data) == 0:
        st.write("No measurement data was uploaded.")
        return

    if st.button(f"Create Dataset from {len(mod_data):,} Abundance Tables"):
        mdata = MuData({
            mod_name: AnnData(mod_df.reindex(index=samples), obs=obs_df)
            for mod_name, mod_df in mod_data
        })
        mdata.pull_obs()

        # Add it to the session state
        set_mdata(mdata, id=id, full=False)
        st.rerun()
