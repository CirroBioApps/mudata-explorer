from collections import defaultdict
import anndata as ad
import hashlib
import json
import muon as mu
import numpy as np
import pandas as pd
import requests
import streamlit as st
from streamlit.errors import StreamlitAPIException
from streamlit.delta_generator import DeltaGenerator
from tempfile import NamedTemporaryFile
from typing import Any, List, Optional, Tuple, Union, Dict
from mudata_explorer.helpers.assets import get_view_by_type, all_view_types # noqa
from mudata_explorer.helpers.assets import get_process_by_type, all_process_types # noqa
from mudata_explorer.helpers.join_kws import join_kws
from mudata_explorer.helpers import mudata, plotting
from mudata_explorer.helpers.save_load import load_history
from mudata_explorer.base.slice import MuDataSlice
from plotly import io


def get_edit_views_flag() -> bool:
    return st.session_state.get("_edit_views_flag", False)


def set_edit_views_flag(val: bool):
    assert isinstance(val, bool)
    st.session_state["_edit_views_flag"] = val


def sidebar_page_links(page_links):
    for path, label, icon in page_links:
        st.sidebar.page_link(
            f"pages/{path}.py",
            label=label,
            icon=icon
        )


def sidebar_edit_views():

    st.sidebar.checkbox(
        "Edit Figures",
        value=get_edit_views_flag(),
        help="Display a set of menus to modify the figures.",
        on_change=update_edit_views,
        key="sidebar_edit_views"
    )


def update_edit_views():
    flag = st.session_state["sidebar_edit_views"]
    if flag != get_edit_views_flag():
        set_edit_views_flag(flag)


def sidebar_load_history():
    if get_mdata() is None:
        return

    if not has_history(exclude=['add_data', 'add_view']):
        return

    if st.sidebar.button("Rerun Analysis"):
        load_history()


@st.cache_resource
def _load_url(url: str):
    res = requests.get(url)
    with NamedTemporaryFile(suffix=".h5mu", delete=True) as tmp:
        with open(tmp.file.name, "wb") as f:
            f.write(res.content)
        return mu.read_h5mu(tmp.file.name)


def load_url(url: str):
    mdata = _load_url(url)
    hydrate_uns(mdata)
    set_mdata(mdata.copy())


def setup_sidebar(
    edit_views=False,
    load_history=False,
    page_layout="centered"
):
    """
    Setup the sidebar with links to all of the pages.
    If edit_views is True, add a checkbox to allow the user to edit the views.
    """
    try:
        st.set_page_config("MuData Explorer", layout=page_layout)
    except StreamlitAPIException:
        st.rerun()

    # If a file link is in the query params
    if st.query_params.get("file"):
        load_url(st.query_params["file"])
        del st.query_params["file"]

    sidebar_page_links([
        ("save_load", "Save / Load", ":material/save:"),
        ("tables", "Tables", ":material/table:"),
        ("processes", "Analysis", ":material/function:"),
        ("views", "Figures", ":material/insert_chart:"),
        ("history", "History", ":material/history:"),
        ("about", "About", ":material/info:")
    ])
    if edit_views:
        sidebar_edit_views()
    if load_history:
        sidebar_load_history()

    plotting.plot_mdata(st.sidebar)


def landing_shortcuts():

    show_shortcuts([
        ("tables", "Upload Tables (*.csv)", ":material/table:"),
        ("about", "About", ":material/info:")
    ])


def show_shortcuts(
    shortcuts: List[tuple],
    ncol=3,
    container: Union[None, DeltaGenerator] = None
):
    cols = None

    for ix, (path, label, icon) in enumerate(shortcuts):
        if ix % ncol == 0:
            cols = (
                st.columns(ncol)
                if container is None
                else container.columns(ncol)
            )

        cols[ix % 3].page_link(
            f"pages/{path}.py",
            label=label,
            icon=icon,
            use_container_width=True
        )


def mdata_to_binary(mdata: mu.MuData) -> bytes:
    # Write out the MuData object to a temporary file
    with NamedTemporaryFile(suffix=".h5mu", delete=True) as tmp:

        # Write the MuData object to the file
        write_h5mu(mdata, tmp.file.name)

        # Get the file object in bytes
        return tmp.read()


def write_h5mu(mdata: mu.MuData, file_name: str):
    # Convert any .uns objects to strings
    dehydrate_uns(mdata)

    # Format as h5mu
    mdata.write(file_name)

    # Convert back to objects
    hydrate_uns(mdata)


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


def jsonify(dat):
    if isinstance(dat, (list, np.ndarray)):
        return [jsonify(val) for val in dat]
    elif isinstance(dat, dict):
        return {kw: jsonify(val) for kw, val in dat.items()}
    else:
        return dat


def validate_json(dat):
    """Validate that an object can be serialized to JSON"""
    try:
        return json.loads(
            json.dumps(
                jsonify(dat),
                sort_keys=True
            )
        )
    except Exception as e:
        print(dat)
        raise ValueError(f"Could not serialize object to JSON: {e}")


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


def duplicate_view(ix: int):
    views = get_views()
    views.insert(ix, views[ix])
    set_views(views)


def add_view(
    view_type: str,
    mdata: Optional[mu.MuData] = None,
    params: Optional[dict] = None
):
    if mdata is None:
        if get_mdata() is None:
            setup_mdata()
    views = get_views(mdata)
    views.append(
        get_view_by_type(view_type).template()
    )
    if params is not None:
        views[-1]["params"] = params
    set_views(views, mdata)


def get_views(mdata: Optional[mu.MuData] = None) -> List[dict]:
    if mdata is None:
        mdata = get_mdata()
    if mdata is None:
        return []
    assert isinstance(mdata, mu.MuData), type(mdata)
    return _json_safe(mdata.uns.get("mudata-explorer-views", []))


def set_views(views, mdata: Optional[mu.MuData] = None):
    use_global = mdata is None
    if use_global:
        mdata = get_mdata()

    # Make sure that the data is JSON serializable
    views = validate_json(views)

    mdata.uns["mudata-explorer-views"] = views
    if use_global:
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


def get_history(exclude=[]) -> List[dict]:
    if not has_mdata():
        return []
    else:
        return [
            h for h in _json_safe(
                get_mdata().uns.get("mudata-explorer-history", [])
            )
            if h.get("process") not in exclude
        ]


def has_history(exclude=[]) -> bool:
    return len(get_history(exclude=exclude)) > 0


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


def set_mdata_hash(mdata_hash: str):
    """Record the hash of the data in the history."""
    add_history(dict(
        process="data_hash",
        params=dict(
            hash=mdata_hash
        ),
        timestamp=get_timestamp()
    ))


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


def _is_axis(axis: int):
    assert axis in [0, 1], f"Unexpected axis: {axis}"


def tree_tables(axis: int) -> List[str]:
    """
    Return a list of all tables in the MuData object.
    axis = 0 for observations, 1 for variables
    """
    _is_axis(axis)

    tables = []

    if has_mdata():

        # If the orientation is to observations,
        # and there is observation metadata
        if axis == 0:
            if get_mdata().obs.shape[1] > 0:
                tables.append("Observation Metadata")

        # Add tables for each modality
        tables.extend([
            join_kws(modality, table)
            for modality in list_modalities()
            for table in list_tables(modality, axis)
        ])

    return tables


def list_tables(modality: str, axis: int):
    if not has_mdata():
        return []
    _is_axis(axis)
    mdata = get_mdata()
    adata: ad.AnnData = mdata.mod[modality]
    tables = ["data"]
    if axis == 0:
        for attr in ["obsm", "obsp"]:
            for slot in getattr(adata, attr).keys():
                tables.append(f"{attr}.{slot}")
    else:
        if adata.var.shape[1] > 0:
            tables.append("metadata")
        for attr in ["varm", "varp"]:
            for slot in getattr(adata, attr).keys():
                tables.append(f"{attr}.{slot}")
    return tables


def list_cnames(tables: Union[str, List[str]], axis=0):

    _is_axis(axis)

    if not has_mdata():
        return []

    mdata = get_mdata()

    if isinstance(tables, str) and len(tables) > 0:
        return get_table_cnames(tables, axis, mdata)

    cnames = []

    for table in tables:
        cnames.extend(
            get_table_cnames(table, axis, mdata)
        )

    # Deduplicate the list of names
    cnames = list(set(cnames))
    cnames.sort()

    return cnames


def get_table_cnames(
    table: str,
    axis: int,
    mdata: mu.MuData
) -> List[str]:

    if table == "Observation Metadata":
        assert axis == 0
        return list(mdata.obs.columns)

    # Otherwise, parse the modality from the name
    assert "." in table, table
    modality, table = table.split(".", 1)
    adata: ad.AnnData = mdata.mod[modality]

    if table == 'metadata':
        cnames = (
            (
                # Note that observations metadata is on mdata.obs
                # while variable metadata is on adata.var
                mdata.obs
                if axis == 0
                else adata.var
            )
            .columns
        )

    elif table == 'data':
        cnames = getattr(
            adata.to_df(),
            "columns" if axis == 0 else "index"
        )

    else:
        prefix, name = table.split(".", 1)
        assert hasattr(adata, prefix), f"Invalid table: {table}"
        cnames = getattr(adata, prefix)[name].columns

    return list(cnames)


def get_dataframe_table(
    modality: str,
    table: str,
    axis: int,
    mdata: Optional[mu.MuData] = None
) -> pd.DataFrame:

    _is_axis(axis)

    # Get the complete set of data
    if mdata is None:
        mdata = get_mdata()

    # Special case for observation metadata
    if axis == 0 and table == "metadata":
        return mdata.obs

    if modality not in mdata.mod:
        return None

    # Get the modality
    adata: ad.AnnData = mdata.mod[modality]

    # Get the table
    if table == "metadata":
        table = mdata.obs if axis == 0 else adata.var
    elif table == "data":
        table = adata.to_df()
        if axis == 1:
            table = table.T
    else:
        assert "." in table, f"Invalid table: {table}"
        attr, kw = table.split(".", 1)
        table = getattr(adata, attr)[kw]

    return table


def join_dataframe_tables(
    tables: List[str],
    axis: int,
    mdata: Optional[mu.MuData] = None
) -> pd.DataFrame:
    """
    Tables from the same modality:
        Joined along the index (if the orientation is to observations)
        Joined along the columns (if the orientation is to variables)
    Tables from different modalities:
        Joined along the columns (if the orientation is to observations)
        Joined along the index (if the orientation is to variables)

    """
    _is_axis(axis)

    # Keep track of which tables are from the same modality
    modality_tables = defaultdict(list)

    if not isinstance(tables, list):
        tables = [tables]

    for table in tables:
        if table == "Observation Metadata":
            assert axis == 0
            modality = 'None'
            df = get_dataframe_table(None, "metadata", axis, mdata=mdata)
        elif '.' in table:
            modality, attr = table.split(".", 1)
            df = get_dataframe_table(modality, attr, axis, mdata=mdata)
        else:
            raise ValueError(f"Could not find table: {table}")

        assert df is not None, f"Could not find table: {table}"
        modality_tables[modality].append(df)

    # Join within each modality
    modalities = {
        modality: (
            pd.concat(tables, axis=axis)
            if len(tables) > 1 else tables[0]
        )
        for modality, tables in modality_tables.items()
    }

    # Join across different modalities
    if len(modalities) > 1:
        df = pd.concat(
            modalities.values(),
            axis=not axis
        )
    elif len(modalities) == 1:
        df: pd.DataFrame = list(modalities.values())[0]
    else:
        return None

    # Drop any rows or columns which are missing in their entirety
    df = df.dropna(axis=0, how="all").dropna(axis=1, how="all")

    return df


def get_dataframe_column(
    mdata: Optional[mu.MuData],
    axis: int,
    table: List[str],
    cname: str
):
    if table is None:
        return

    _is_axis(axis)

    df = join_dataframe_tables(
        table,
        axis,
        mdata=mdata
    )
    # Just get the column of interest
    df = df[cname]

    # If there are multiple columns with the same name
    if len(df.shape) > 1:
        df = df.iloc[:, 0]

    return df


def get_dat_hash(mdata: Optional[mu.MuData] = None):
    if mdata is None and has_mdata() is False:
        return None, None, None

    if mdata is None:
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
        with container.expander(
            f"**Provenance: '{loc.address}'**"
        ):
            st.write(f"**Analysis:** {prov['process']}")
            st.write(f"**Timestamp:** {prov['timestamp']}")
            st.write("**SDK Configuration:**")
            st.code(process_sdk_snippet(prov))

            if isinstance(prov.get("figures"), list):
                if len(prov.get("figures")) > 0:
                    st.write("Supporting Figures")
                for fig_json in prov["figures"]:
                    fig = io.from_json(fig_json)
                    st.plotly_chart(fig)
        return True
    return False


def nest_params(params: dict):
    output = dict()
    for key, value in params.items():
        if '.' in key:
            keys = key.split('.', 1)
            if keys[0] not in output:
                output[keys[0]] = dict()
            output[keys[0]][keys[1]] = value
        else:
            output[key] = value
    return {
        key: nest_params(value) if isinstance(value, dict) else value
        for key, value in output.items()
    }


def process_sdk_snippet(prov: dict):
    assert "process" in prov.keys()
    assert "params" in prov.keys()
    params = nest_params(prov["params"])
    params_str = (
        json.dumps(params, indent=4)
        .replace('false', 'False')
        .replace('true', 'True')
        .replace('null', 'None')
        .replace("\n", "\n    ")
    )
    return f"""process.{prov['process'].replace('-', '_')}(
    mdata,
    **{params_str}
)
"""


def get_supp_figs() -> List[Tuple[str, dict]]:
    """Return the list of figures which are stored in the provenance."""

    return [
        f"{loc}:{ix}"
        for loc, prov in get_provenance().items()
        if prov.get("figures") is not None
        for ix, fig in enumerate(prov["figures"])
        if fig is not None
    ]
