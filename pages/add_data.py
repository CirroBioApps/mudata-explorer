from copy import copy
from mudata_explorer import app

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


def add_mudata(obs: pd.DataFrame, mod: pd.DataFrame, mod_name: str):
    overlap = set(obs.index).intersection(set(mod.index))
    obs = obs.loc[list(overlap)]
    mod = mod.loc[list(overlap)]

    adata = ad.AnnData(obs=obs, X=mod)

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


if __name__ == "__main__":

    app.setup_pages()

    container = st.container()

    if app.get_mdata() is not None:
        container.page_link("pages/summarize.py", label="**View Existing Data**")

    container.write("#### Add Data")

    # Get input from the user
    obs_file = container.file_uploader(
        "Observation Metadata (.obs)",
        help="Provide a CSV/TSV where the first column is a unique identifier for each observation.", # noqa
    )
    obs = read_table(obs_file, container)

    mod_name = container.text_input(
        "Measurement Name",
        value="Measurement",
        help="Enter the name of the measurement.",
    )

    mod_file = container.file_uploader(
        "Measurement Data (.X)",
        help="Provide a CSV/TSV where the first column is a unique identifier for each observation.", # noqa
    )
    mod = read_table(mod_file, container)

    if obs is None:
        container.write("No metadata uploaded")

    else:
        obs = sanitize_types(obs, container, keep_str=True)
        container.write(f"Metadata: {obs.shape[0]:,} rows x {obs.shape[1]:,} columns.") # noqa

    if mod is None:
        container.write("No measurement data uploaded.")

    else:

        mod = sanitize_types(mod, container)
        container.write(f"{mod_name}: {mod.shape[0]:,} rows x {mod.shape[1]:,} columns.") # noqa

    if mod is not None and obs is not None:

        # Find the overlapping index labels between the two tables
        overlap = set(obs.index).intersection(set(mod.index))

        container.write(f"No. of overlapping observations: {len(overlap):,}")

        if len(overlap) > 0:

            # Figure out if all of the columns should be added
            if container.checkbox("Use all columns", value=True):
                columns = mod.columns.values
            else:
                columns = container.multiselect(
                    "Select columns",
                    mod.columns.values,
                    default=mod.columns.values
                )

            if container.button(
                "Add to MuData",
                on_click=add_mudata,
                args=(
                    obs,
                    mod.reindex(columns=columns),
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
