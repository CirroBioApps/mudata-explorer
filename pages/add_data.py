from copy import copy
from typing import Union
from mudata_explorer import app
from mudata_explorer.helpers import plotting

import anndata as ad
import muon as mu
import pandas as pd
from streamlit.delta_generator import DeltaGenerator
import streamlit as st


def sanitize_types(
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


def read_table(file, container: DeltaGenerator):

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

    # Set the first column as the index
    df = df.set_index(df.columns[0])

    # Return the data
    return df


def add_mudata(adata: ad.AnnData, mod_name: str):

    mdata = app.get_mdata()
    if mdata is None:
        mdata = mu.MuData({
            mod_name: adata
        })

    else:

        uns = copy(mdata.uns)
        mdata = mu.MuData({
            **{
                kw: adata
                for kw, adata in mdata.mod.items()
                if kw != '_blank'
            },
            **{
                mod_name: adata
            }
        })
        for kw, val in uns.items():
            mdata.uns[kw] = val

    mdata.update_obs()
    app.set_mdata(mdata)

    event = dict(
        timestamp=app.get_timestamp(),
        process="add_data",
        params=dict(name=mod_name),
        updated_keys=["X"]
    )
    app.add_provenance(
        mod_name,
        "X",
        None,
        event
    )
    app.add_history(event)


def _get_obs(container: DeltaGenerator):

    # If there is no data, ask the user for observation metadata
    if app.get_mdata() is None:

        # Get a metadata table from the user
        return _get_table(
            container,
            "Observation Metadata (.obs)",
            help="Provide a CSV/TSV where the first column is a unique identifier for each observation.", # noqa
            keep_str=True
        )

    else:
        return None


def _get_name(container: DeltaGenerator):

    mod_name = container.text_input(
        "Measurement Name",
        value="Measurement",
        help="Enter the name of the measurement.",
    )
    if "." in mod_name:
        container.error("The measurement name cannot contain a period.")
        return
    else:
        if mod_name in app.list_modalities():
            container.error("A measurement with this name already exists.")
            return
        return mod_name


def _get_df(container: DeltaGenerator):
    return _get_table(
        container,
        "Measurement Data (.X)",
        help="Provide a CSV/TSV where the first column is a unique identifier for each observation.", # noqa
    )


def _get_var(container: DeltaGenerator):

    return _get_table(
        container,
        "Variable Metadata (.var)",
        help="Provide a CSV/TSV where the first column is a unique identifier for each variable.", # noqa
        keep_str=True
    )


def _get_table(
    container: DeltaGenerator,
    label: str,
    help=None,
    keep_str=False
):

    # Get input from the user
    file = container.file_uploader(
        label,
        help=help
    )
    df = read_table(file, container)
    if df is not None:
        df = sanitize_types(df, container, keep_str=keep_str)
    return df


def _df_desc_str(name: str, df: pd.DataFrame):
    return f"{name}: {df.shape[0]:,} rows x {df.shape[1]:,} columns."


def _merge_data(
    obs: Union[pd.DataFrame, None],
    df: pd.DataFrame,
    var: Union[pd.DataFrame, None],
    container: DeltaGenerator
):

    if df is None:
        container.write("No measurement data uploaded.")
        return None

    if obs is not None:
        container.write(_df_desc_str("Observation Metadata", obs))

    container.write(_df_desc_str("Measurement Data", df))

    if var is not None:
        container.write(_df_desc_str("Variable Metadata", var))

    # Report if there are any rows in obs or var which are not in the data
    missing_str = "There are {} in the metadata which are not in the data."

    if obs is not None:
        # Find if there are any values in obs.index which are not in df.index
        missing = set(obs.index).difference(set(df.index))
        if len(missing) > 0:
            container.write(missing_str.format(f"{len(missing):,} observations"))

        # Only keep the observations which are in the data
        obs = obs.reindex(index=df.index)

    if var is not None:
        # Find if there are any values in var.index which are not in df.columns
        missing = set(var.index).difference(set(df.columns))
        if len(missing) > 0:
            container.write(missing_str.format(f"{len(missing):,} variables"))

        # Only keep the variables which are in the data
        var = var.reindex(index=df.columns)

    # Figure out if all of the columns should be added
    if not container.checkbox("Use all columns (variables)", value=True):
        df = df.reindex(
            columns=container.multiselect(
                "Select columns (variables)",
                df.columns.values,
                default=df.columns.values
            )
        )

    # Figure out if all of the rows should be added
    if not container.checkbox("Use all rows (observations)", value=True):
        df = df.reindex(
            index=container.multiselect(
                "Select rows (observations)",
                df.index.values,
                default=df.index.values
            )
        )

    # Make an AnnData object
    return ad.AnnData(
        X=df,
        obs=obs,
        var=var
    )


if __name__ == "__main__":

    app.setup_pages()

    container = st.container()

    if app.get_mdata() is not None:
        container.page_link("pages/summarize.py", label="**View Existing Data**")

    container.write("#### Add Data")

    # Show the dataset
    plotting.plot_mdata(container)

    # Get the observation metadata (if needed)
    obs = _get_obs(container)

    # Get the measurement data
    mod_name = _get_name(container)
    df = _get_df(container)

    # Get the variable annotations (if any)
    var = _get_var(container)

    # Merge into a single AnnData object
    adata = _merge_data(obs, df, var, container)

    if container.button(
        "Add to MuData",
        on_click=add_mudata,
        args=(
            adata,
            mod_name
        )
    ):
        container.write("Data added to MuData.")
        app.show_shortcuts(
            [
                ("summarize", ":book: Inspect Uploaded Data"),
                ("views", ":bar_chart: View Data"),
                ("processes", ":running: Run Processes")
            ],
            container=container
        )
