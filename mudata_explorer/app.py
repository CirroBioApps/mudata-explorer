import anndata as ad
import hashlib
import json
import muon as mu
import numpy as np
import pandas as pd
import streamlit as st
from streamlit.delta_generator import DeltaGenerator
from tempfile import NamedTemporaryFile
from typing import Any, List, Union, Dict
from mudata_explorer.helpers import get_view_by_type
from mudata_explorer.helpers.join_kws import join_kws


def page_links():
    return [
        ("summarize", "Summarize"),
        ("views", "Views"),
        ("processes", "Process Data"),
        ("add_data", "Add Data"),
        ("save_load", "Save/Load"),
        ("history", "History"),
        ("settings", "Settings"),
        ("about", "About")
    ]


def setup_pages():
    st.set_page_config("MuData Explorer", layout="centered")
    st.sidebar.title("MuData Explorer")

    for path, label in page_links():
        st.sidebar.page_link(
            f"pages/{path}.py",
            label=label
        )


def landing_shortcuts():

    show_shortcuts([
        ("add_data", ":page_facing_up: Upload Tables (*.csv)"),
        ("save_load", ":open_file_folder: Load Dataset (*.h5mu)"),
        ("about", ":information_source: About")
    ])


def show_shortcuts(
    shortcuts: List[tuple],
    ncol=3,
    container: Union[None, DeltaGenerator] = None
):
    cols = None

    for ix, (path, label) in enumerate(shortcuts):
        if ix % ncol == 0:
            cols = (
                st.columns(ncol)
                if container is None
                else container.columns(ncol)
            )

        cols[ix % 3].page_link(
            f"pages/{path}.py",
            label=label,
            use_container_width=True
        )


def mdata_to_binary(mdata: mu.MuData) -> bytes:
    # Write out the MuData object to a temporary file
    with NamedTemporaryFile(suffix=".h5mu", delete=True) as tmp:

        # Convert any .uns objects to strings
        dehydrate_uns(mdata)

        # Format as h5mu
        mdata.write(tmp.file.name)

        # Convert back to objects
        hydrate_uns(mdata)

        # Get the file object in bytes
        return tmp.read()


def read_h5mu(h5mu_file):

    with NamedTemporaryFile(suffix=".h5mu", delete=True) as tmp:
        with open(tmp.file.name, "wb") as f:
            f.write(h5mu_file.getvalue())
        mdata = mu.read_h5mu(tmp.file.name)

    hydrate_uns(mdata)

    return mdata


def hash_dat(dat, n: Union[int, None] = 16):
    """Compute the hash of an object."""
    hash = hashlib.sha256()
    hash.update(dat)
    hex = hash.hexdigest()
    if n is not None:
        hex = hex[:n]
    return hex


def get_mdata() -> Union[None, mu.MuData]:
    mdata = st.session_state.get("mdata", None)
    if mdata is not None:
        assert isinstance(mdata, mu.MuData)
        if "mudata-explorer-process" not in mdata.uns.keys():
            mdata.uns["mudata-explorer-process"] = {
                "category": None,
                "type": None,
                "params": {}
            }
    return mdata


def set_mdata(
    mdata: mu.MuData,
    timestamp: Union[None, str] = None,
    process: Union[None, str] = None,
    params: Union[None, Dict[str, Any]] = None,
    updated_keys: Union[List[str], str] = None
):
    """
    Set the MuData object

    Optionally include process information, which will be saved
    as history and provenance information.

    Optional Args:

        timestamp = str e.g. app.get_timestamp()
        process = str
        params = dict
        updated_keys = List[str] e.g. ["rna.X", "rna.obs", "rna.var"]
    """

    assert isinstance(mdata, mu.MuData)
    if (
        timestamp is not None or
        process is not None or
        params is not None
    ):
        event = dict(
            timestamp=timestamp,
            process=process,
            params=params,
            updated_keys=updated_keys
        )

        # Add the event to the history
        add_history(event, mdata)

        # If any updated keys were provided
        if isinstance(updated_keys, str):
            updated_keys = [updated_keys]
        if isinstance(updated_keys, list):
            for kw in updated_keys:
                add_provenance(kw, event, mdata)

    st.session_state["mdata"] = mdata


def setup_mdata():
    mdata = mu.MuData({
        '_blank': ad.AnnData(
            X=np.array([[]]),
            obs=[],
            var=[]
        )
    })
    mdata.uns["mudata-explorer-views"] = []
    mdata.uns["mudata-explorer-process"] = {}
    mdata.uns["mudata-explorer-settings"] = {}
    mdata.uns["mudata-explorer-history"] = []
    mdata.uns["mudata-explorer-provenance"] = {}
    set_mdata(mdata)


def validate_json(dat):
    """Validate that an object can be serialized to JSON"""
    return json.loads(json.dumps(dat, sort_keys=True))


def get_process() -> dict:
    mdata = get_mdata()
    if mdata is None:
        return {}
    assert isinstance(mdata, mu.MuData), type(mdata)
    return mdata.uns.get("mudata-explorer-process", {})


def set_process(process: dict) -> None:
    mdata = get_mdata()
    assert mdata is not None
    assert isinstance(mdata, mu.MuData), type(mdata)
    mdata.uns["mudata-explorer-process"] = process
    set_mdata(mdata)


def update_process_on_change(kw) -> None:
    val = st.session_state[f"process-{kw}"]
    update_process(kw, val)


def update_process(kw, val) -> None:
    process = get_process()
    process[kw] = val
    set_process(process)


def delete_view(ix: int):
    views = get_views()
    views.pop(ix)
    set_views(views)


def add_view(view_type: str):
    if get_mdata() is None:
        setup_mdata()
    views = get_views()
    views.append(
        get_view_by_type(view_type).template()
    )
    set_views(views)


def get_views() -> List[dict]:
    mdata = get_mdata()
    if mdata is None:
        return []
    assert isinstance(mdata, mu.MuData), type(mdata)
    return mdata.uns.get("mudata-explorer-views", [])


def set_views(views):
    mdata = get_mdata()

    # Make sure that the data is JSON serializable
    views = validate_json(views)

    mdata.uns["mudata-explorer-views"] = views
    set_mdata(mdata)


def get_settings() -> dict:
    mdata = get_mdata()
    if mdata is None:
        settings = {}
    else:
        settings = get_mdata().uns.get("mudata-explorer-settings", {})

    for kw, val in dict(
        editable=True
    ).items():
        if kw not in settings:
            settings[kw] = val
    return settings


def set_settings(settings: dict):
    mdata = get_mdata()
    assert isinstance(mdata, mu.MuData), type(mdata)

    # Make sure that the data is JSON serializable
    settings = validate_json(settings)

    if json.dumps(settings) != json.dumps(get_settings()):

        mdata.uns["mudata-explorer-settings"] = settings
        set_mdata(mdata)
        st.experimental_rerun()


def get_history() -> List[dict]:
    mdata = get_mdata()

    if mdata is None:
        return []

    assert isinstance(mdata, mu.MuData), type(mdata)

    return mdata.uns.get("mudata-explorer-history", [])


def set_history(history: dict):
    mdata = get_mdata()
    assert isinstance(mdata, mu.MuData)

    # Make sure that the data is JSON serializable
    history = validate_json(history)

    mdata.uns["mudata-explorer-history"] = history
    set_mdata(mdata)


def add_history(event: dict):
    mdata = get_mdata()
    if mdata is None:
        return
    history = get_history()
    history.insert(0, event)
    set_history(history)


def get_provenance() -> Dict[str, dict]:
    mdata = get_mdata()

    if mdata is None:
        return {}
    assert isinstance(mdata, mu.MuData), type(mdata)

    return mdata.uns.get("mudata-explorer-provenance", {})


def query_provenance(
    mod_name: str,
    slot: str,
    kw: str
) -> Union[None, dict]:
    key = join_kws(mod_name, slot, kw)
    provenance = get_provenance()
    return provenance.get(key, None)


def set_provenance(provenance: dict):
    mdata = get_mdata()
    assert isinstance(mdata, mu.MuData), type(mdata)

    # Make sure that the data is JSON serializable
    provenance = validate_json(provenance)

    mdata.uns["mudata-explorer-provenance"] = provenance
    set_mdata(mdata)


def add_provenance(
    mod_name: str,
    slot: str,
    kw: Union[str, None],
    event: dict
):
    provenance = get_provenance()
    provenance_key = join_kws(mod_name, slot, kw)
    provenance[provenance_key] = event
    set_provenance(provenance)


def hydrate_uns(mdata: mu.MuData):
    prefix = "mudata-explorer-"
    for suffix in ["views", "history", "provenance", "settings"]:
        kw = f"{prefix}{suffix}"
        if kw in mdata.uns:
            if isinstance(mdata.uns[kw], str):
                mdata.uns[kw] = json.loads(mdata.uns[kw])


def dehydrate_uns(mdata: mu.MuData):
    prefix = "mudata-explorer-"
    for suffix in ["views", "history", "provenance", "settings"]:
        kw = f"{prefix}{suffix}"
        if kw in mdata.uns:
            if not isinstance(mdata.uns[kw], str):
                mdata.uns[kw] = json.dumps(mdata.uns[kw], sort_keys=True)


def make_modality_df(mdata: mu.MuData, mod_name: str) -> pd.DataFrame:
    """
    Make a single DataFrame which combines information from the modality.
    If there are any conflicting column names, they will be suffixed with
    the slot they came from.
    """
    adata: ad.AnnData = mdata.mod[mod_name]

    # Collect all of the DataFrames
    all_dfs = {}

    # Data from the X slot
    all_dfs["data"] = adata.to_df()

    # Metadata
    all_dfs["metadata"] = adata.obs

    # For each obsm slot
    for kw, slot in adata.obsm.items():
        all_dfs[kw] = pd.DataFrame(slot, index=adata.obs.index)

    # Find each column name which is repeated across the DataFrames
    all_cols = []
    for df in all_dfs.values():
        all_cols.extend(df.columns.values)
    all_cols = pd.Series(all_cols)
    dup_cols = all_cols[all_cols.duplicated()].unique()

    # For each duplicated column name, add the slot name as a suffix
    for col in dup_cols:
        for slot, df in all_dfs.items():
            if col in df.columns:
                df.rename(columns={col: f"{col}_{slot}"}, inplace=True)

    # Concatenate all of the DataFrames
    df = pd.concat(all_dfs.values(), axis=1)

    return df


def get_timestamp():
    return str(pd.Timestamp.now())


def list_modalities():
    mdata = get_mdata()
    if mdata is None:
        return []
    return list(mdata.mod.keys())


def tree_tables(orientation) -> List[str]:
    """Return a list of all tables in the MuData object."""
    return [
        join_kws(modality, table)
        for modality in list_modalities()
        for table in list_tables(modality, orientation)
    ]


def list_tables(modality: str, orientation: str):
    assert orientation in ["observations", "variables"], orientation
    mdata = get_mdata()
    if mdata is None:
        return []
    adata: ad.AnnData = mdata.mod[modality]
    tables = ["metadata", "data"]
    if orientation == "observations":
        for attr in ["obsm", "obsp"]:
            for slot in getattr(adata, attr).keys():
                tables.append(f"{attr}.{slot}")
    else:
        for attr in ["varm", "varp"]:
            for slot in getattr(adata, attr).keys():
                tables.append(f"{attr}.{slot}")
    return tables


def list_cnames(modality: str, table: str, orientation="observations"):
    if orientation == "observations":
        attr = "columns"
    else:
        assert orientation == "variables"
        attr = "index"

    mdata = get_mdata()
    if mdata is None:
        return []
    adata: ad.AnnData = mdata.mod[modality]
    if table == 'metadata':
        return list(getattr(adata.obs, attr))
    elif table == 'data':
        return list(getattr(adata.to_df(), attr))
    else:
        prefix, name = table.split(".", 1)
        assert hasattr(adata, prefix)
        return list(getattr(adata, prefix)[name].columns)
    raise ValueError(f"Invalid table: {table}")


def get_dataframe_table(modality: str, table: str) -> pd.DataFrame:

    # Get the complete set of data
    mdata = get_mdata()

    if modality not in mdata.mod:
        return None

    # Get the modality
    adata: ad.AnnData = mdata.mod[modality]

    # Get the table
    if table == "metadata":
        table = adata.obs
    elif table == "data":
        table = adata.to_df()
    elif table.startswith("obsm."):
        table = adata.obsm[table[5:]]
    elif table.startswith("obsp."):
        table = adata.obsp[table[5:]]
    else:
        raise ValueError(f"Invalid table: {table}")

    return table


def get_dataframe_column(
    modality: str,
    table: str,
    cname: str
):

    table = get_dataframe_table(modality, table)

    # Get the column
    if cname not in table.columns:
        return None

    return table[cname]


def get_dat_hash():
    mdata = get_mdata()
    if mdata is None:
        return None, None

    # Convert the MuData object to binary
    dat = mdata_to_binary(mdata)

    # Compute the hash of the file
    hash = hash_dat(dat)

    return dat, hash


def save_annot(
    mdata: mu.MuData,
    modality: str,
    slot: str,
    dest_key: str,
    column_dat: Union[pd.Series, pd.DataFrame],
    params: dict,
    process_type: str
):

    # Add the PCA coordinates to the obs/var slot
    getattr(mdata.mod[modality], slot)[dest_key] = column_dat

    # Make a record of the process
    event = dict(
        process=process_type,
        params=params,
        timestamp=get_timestamp(),
        updated_keys=dest_key
    )

    # Save the results

    # Update the MuData object
    set_mdata(mdata)

    # Add it to the history
    add_history(event)

    # Mark the source of the table which was added
    add_provenance(
        modality,
        slot,
        dest_key,
        event
    )


def show_provenance(mdata, modality, slot, dest_key, container):

    # If the key already exists
    if dest_key in getattr(mdata.mod[modality], slot).keys():
        prov = query_provenance(modality, slot, dest_key)
        if prov is not None:
            container.write(f"**Data currently in '{slot}/{dest_key}'**")
            container.write(prov)
            container.write(
                f"> 'Run' will overwrite existing data in '{slot}/{dest_key}'."
            )
