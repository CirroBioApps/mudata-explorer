import anndata as ad
import hashlib
import json
import muon as mu
import pandas as pd
import streamlit as st
from streamlit.delta_generator import DeltaGenerator
from tempfile import NamedTemporaryFile
from typing import Any, List, Optional, Union, Dict
from mudata_explorer.helpers import get_view_by_type
from mudata_explorer.helpers.join_kws import join_kws
from mudata_explorer.helpers import mudata, plotting, save_load
from mudata_explorer.base.slice import MuDataSlice
from plotly import io


def sidebar_page_links(page_links):
    for path, label in page_links:
        st.sidebar.page_link(
            f"pages/{path}.py",
            label=label
        )


def sidebar_edit_views():
    if get_mdata() is None:
        return

    settings = get_settings()

    settings["editable"] = st.sidebar.checkbox(
        "Edit Figures",
        value=settings.get("editable", True),
        help="Display a set of menus to modify the figures.",
        on_change=update_edit_views,
        key="sidebar_edit_views"
    )


def update_edit_views():
    flag = st.session_state["sidebar_edit_views"]
    settings = get_settings()
    if flag != settings["editable"]:
        settings["editable"] = flag
        set_settings(settings)


def setup_sidebar(edit_views=False):
    """
    Setup the sidebar with links to all of the pages.
    If edit_views is True, add a checkbox to allow the user to edit the views.
    """
    st.set_page_config("MuData Explorer", layout="centered")
    # st.sidebar.title("MuData Explorer")

    sidebar_page_links([
        ("tables", "Tables"),
        ("processes", "Analysis"),
        ("views", "Figures"),
        ("history", "History"),
        # ("settings", "Settings"),
        ("about", "About")
    ])
    if edit_views:
        sidebar_edit_views()

    load_cont, save_cont = st.sidebar.columns(2)
    if has_mdata():
        if save_cont.button("Save File", use_container_width=True):
            save_load.download_button()
    if load_cont.button("Load File", use_container_width=True):
        save_load.upload_button()

    plotting.plot_mdata(st.sidebar)


def landing_shortcuts():

    show_shortcuts([
        ("tables", ":page_facing_up: Upload Tables (*.csv)"),
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
    if "mdata" not in st.session_state:
        mdata = None
    else:
        mdata = st.session_state["mdata"]

    if mdata is not None:
        assert isinstance(mdata, mu.MuData)
        if "mudata-explorer-process" not in mdata.uns.keys():
            mdata.uns["mudata-explorer-process"] = {
                "category": None,
                "type": None,
                "params": {}
            }
    return mdata


def has_mdata() -> bool:
    return len(list_modalities()) > 0


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

    assert isinstance(mdata, mu.MuData), type(mdata)
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
    mdata = mudata.setup_mdata()
    set_mdata(mdata)


def validate_json(dat):
    """Validate that an object can be serialized to JSON"""
    return json.loads(json.dumps(dat, sort_keys=True))


def get_process() -> dict:
    if not has_mdata():
        return {}
    mdata = get_mdata()
    assert isinstance(mdata, mu.MuData), type(mdata)
    return _json_safe(mdata.uns.get("mudata-explorer-process", {}))


def set_process(process: dict) -> None:
    if not has_mdata():
        setup_mdata()

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
    if not has_mdata():
        return []
    mdata = get_mdata()
    assert isinstance(mdata, mu.MuData), type(mdata)
    return _json_safe(mdata.uns.get("mudata-explorer-views", []))


def set_views(views):
    mdata = get_mdata()

    # Make sure that the data is JSON serializable
    views = validate_json(views)

    mdata.uns["mudata-explorer-views"] = views
    set_mdata(mdata)


def _json_safe(obj: Union[str, dict]):
    if isinstance(obj, str):
        return json.loads(obj)
    else:
        return obj


def get_settings() -> dict:
    if not has_mdata():
        settings = {}
    else:
        settings = _json_safe(
            get_mdata().uns.get("mudata-explorer-settings", {})
        )

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
    if not has_mdata():
        return []
    else:
        return _json_safe(
            get_mdata().uns.get("mudata-explorer-history", [])
        )


def set_history(history: dict):
    mdata = get_mdata()
    assert isinstance(mdata, mu.MuData)

    # Make sure that the data is JSON serializable
    history = validate_json(history)

    mdata.uns["mudata-explorer-history"] = history
    set_mdata(mdata)


def add_history(event: dict):
    if not has_mdata():
        return
    history = get_history()
    history.insert(0, event)
    set_history(history)


def get_provenance() -> Dict[str, dict]:
    if not has_mdata():
        return {}

    mdata = get_mdata()
    assert isinstance(mdata, mu.MuData), type(mdata)

    return _json_safe(mdata.uns.get("mudata-explorer-provenance", {}))


def query_provenance(loc: MuDataSlice) -> Union[None, dict]:
    provenance = get_provenance()
    return provenance.get(loc.dehydrate(), None)


def set_provenance(provenance: dict):
    mdata = get_mdata()
    assert isinstance(mdata, mu.MuData), type(mdata)

    # Make sure that the data is JSON serializable
    provenance = validate_json(provenance)

    mdata.uns["mudata-explorer-provenance"] = provenance
    set_mdata(mdata)


def add_provenance(
    loc: MuDataSlice,
    event: dict
):
    provenance = get_provenance()
    provenance[loc.dehydrate()] = event
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
    mdata: mu.MuData = get_mdata()
    if mdata is None:
        return []
    mods = [
        mod
        for mod in mdata.mod.keys()
        if not mod.startswith("_")
    ]
    return mods


def tree_tables(orientation) -> List[str]:
    """Return a list of all tables in the MuData object."""
    return [
        join_kws(modality, table)
        for modality in list_modalities()
        for table in list_tables(modality, orientation)
    ]


def list_tables(modality: str, orientation: str):
    if not has_mdata():
        return []
    assert orientation in ["observations", "variables"], orientation
    mdata = get_mdata()
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

    msg = f"Unexpected orientation: {orientation}"
    assert orientation in ["observations", "variables"], msg

    if not has_mdata():
        return []

    mdata = get_mdata()
    adata: ad.AnnData = mdata.mod[modality]

    if table == 'metadata':
        cnames = (
            (
                # Note that observations metadata is on mdata.obs
                # while variable metadata is on adata.var
                mdata.obs
                if orientation == "observations"
                else adata.var
            )
            .columns
        )

    elif table == 'data':
        cnames = getattr(
            adata.to_df(),
            "columns" if orientation == "observations" else "index"
        )

    else:
        prefix, name = table.split(".", 1)
        assert hasattr(adata, prefix), f"Invalid table: {table}"
        cnames = getattr(adata, prefix)[name].columns

    return list(cnames)


def get_dataframe_table(
    modality: str,
    table: str,
    orientation: str
) -> pd.DataFrame:

    msg = f"Unexpected orientation: {orientation}"
    assert orientation in ["observations", "variables"], msg

    # Get the complete set of data
    mdata = get_mdata()

    if modality not in mdata.mod:
        return None

    # Get the modality
    adata: ad.AnnData = mdata.mod[modality]

    # Get the table
    if table == "metadata":
        table = mdata.obs if orientation == "observations" else adata.var
    elif table == "data":
        table = adata.to_df()
        if orientation == "variables":
            table = table.T
    else:
        assert "." in table, f"Invalid table: {table}"
        attr, kw = table.split(".", 1)
        table = getattr(adata, attr)[kw]

    return table


def get_dataframe_column(
    orientation: str,
    modality: str,
    table: str,
    cname: str
):

    # Parse the slot, attr, and subattr from the table name
    # Get the table
    subattr = None
    if table == "metadata":
        if orientation == "observations":
            slot = "obs"
        else:
            slot = "var"
        attr = cname
    elif table == "data":
        slot = "X"
        attr = cname
    elif "." in table:
        slot, attr = table.split(".", 1)
        subattr = cname
    else:
        raise ValueError(f"Invalid table: {table}")

    # Define the slice of the data
    slice = MuDataSlice(
        modality=modality,
        orientation=orientation[:3],
        slot=slot,
        attr=attr,
        subattr=subattr
    )

    # Return the data
    return slice.dataframe(get_mdata())


def get_dat_hash():
    if not has_mdata():
        return None, None, None

    mdata = get_mdata()

    # Convert the MuData object to binary
    dat = mdata_to_binary(mdata)

    # Compute the hash of the file
    hash = hash_dat(dat)

    # Compute the size of the file
    size = len(dat) / 1024

    # Format the size as a string
    if size < 1024:
        size = f"{size:.2f} KB"
    else:
        size = f"{size/1024:.2f} MB"

    return dat, hash, size


def save_annot(
    mdata: mu.MuData,
    loc: MuDataSlice,
    column_dat: Union[pd.Series, pd.DataFrame],
    params: dict,
    process_type: str,
    figures: Optional[List[dict]]
):

    # Write the data to the specified address
    loc.write(mdata, column_dat)

    # Make a record of the process
    event = dict(
        process=process_type,
        params=params,
        timestamp=get_timestamp(),
        loc=loc.params,
        figures=figures
    )

    # Save the results

    # Update the MuData object
    set_mdata(mdata)

    # Add it to the history
    add_history(event)

    # Mark the source of the table which was added
    add_provenance(loc, event)


def show_provenance(loc: MuDataSlice, container: DeltaGenerator):

    prov = query_provenance(loc)
    if prov is not None:
        container.write(
            f"**Data currently in '{loc.address}'**"
        )
        container.write({
            kw: val
            for kw, val in prov.items()
            if kw != "figures"
        })
        if isinstance(prov.get("figures"), list):
            if len(isinstance(prov.get("figures"), list)) > 0:
                container.write("Supporting Figures")
            for fig_json in prov["figures"]:
                fig = io.from_json(fig_json)
                container.plotly_chart(fig)
        container.write(
            f"> 'Run' will overwrite existing data in '{loc.address}'."
        )
