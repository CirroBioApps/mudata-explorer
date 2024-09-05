from copy import copy
from collections import defaultdict
from mudata_explorer.base.slice import MuDataSlice
from mudata_explorer.helpers.io import json_safe, validate_json
from mudata_explorer.helpers.join_kws import join_kws
from mudata_explorer.helpers.timestamp import get_timestamp
from typing import Union, Dict, Any, List, Tuple, Optional
import anndata as ad
import muon as mu
import numpy as np
import pandas as pd
import streamlit as st


"""
The saving and loading of a MuData object in the streamlit
session state involves splitting up each element of the
object so that callback invalidation is limited to the smallest
necessary scope.

The top-level key will be mdata-{id}, where
the default id is 'main'. This makes it possible to store
multiple MuData elements in the session state.

Everything apart from the mdata-* elements
of .uns will be saved to mdata-{id}-data.

The other attributes of .uns will be saved to
mdata-{id}-{key} as appropriate.

The one exception will be views, which will be saved
using:

- mdata-{id}-n_views: int
- mdata-{id}-views-{ix}

A boolean key will be used to note whether any data
exists: mdata-{id}-exists
"""

#####################
# GET / SET METHODS #
#####################

def _session_key(attribute: str, id="main"):
    return f"mdata-{id}-{attribute}"


def _get_mdata_elem(
    attribute: str,
    default=None,
    id="main"
):

    return copy(st.session_state.get(
        _session_key(attribute, id=id),
        default=default
    ))


def _set_mdata_elem(
    attribute: str,
    value=None,
    id="main"
):

    st.session_state[_session_key(attribute, id=id)] = copy(value)


def get_mdata_exists(id="main") -> bool:
    return _get_mdata_elem("exists", default=False, id=id)


def set_mdata_exists(exists: bool, id="main"):
    assert isinstance(exists, bool)
    _set_mdata_elem("exists", value=exists, id=id)


def get_view(ix: int, id="main") -> dict:
    view = _get_mdata_elem(f"views-{ix}", id=id)
    assert view is not None, f"No view exists with ix={ix}"
    return view


def set_view(ix: int, view: dict, id="main"):
    _set_mdata_elem(f"views-{ix}", value=view, id=id)


def get_views(id="main") -> List[dict]:
    return [
        get_view(ix, id=id)
        for ix in range(_get_mdata_elem("n_views", default=0, id=id))
    ]


def set_views(views: List[dict], id="main"):
    assert isinstance(views, list)
    _set_mdata_elem("n_views", value=len(views), id=id)
    for ix, view in enumerate(views):
        assert isinstance(view, dict)
        set_view(ix, view, id=id)


def get_process(id="main") -> dict:
    default_process = {"category": None, "type": None, "params": {}}
    return json_safe(
        _get_mdata_elem("process", default=default_process, id=id)
    )


def set_process(process, id="main"):
    _set_mdata_elem("process", value=process, id=id)


def get_settings(id="main") -> dict:
    return json_safe(_get_mdata_elem("settings", default={}, id=id))


def set_settings(settings, id="main"):
    assert isinstance(settings, dict)
    # Make sure that the data is JSON serializable
    settings = validate_json(settings)
    _set_mdata_elem("settings", value=settings, id=id)


def get_history(id="main", exclude=[]) -> List[dict]:
    return [
        h for h in json_safe(
            _get_mdata_elem("history", default=[], id=id)
        )
        if h.get("process") not in exclude
    ]


def has_history(id="main", exclude=[]) -> bool:
    return len(get_history(id=id, exclude=exclude)) > 0


def set_history(history: list, id="main"):
    assert isinstance(history, list)

    # Make sure that the data is JSON serializable
    history = validate_json(history)

    _set_mdata_elem("history", value=history, id=id)


def add_history(event: dict, id="main"):
    history = get_history(id=id)
    history.insert(0, event)
    set_history(history, id=id)


def get_provenance(id="main") -> Dict[str, dict]:
    return json_safe(
        _get_mdata_elem("provenance", default={}, id=id)
    )

def query_provenance(loc: MuDataSlice, id="main") -> Union[None, dict]:
    provenance = get_provenance(id=id)
    return provenance.get(loc.dehydrate(), None)


def set_provenance(provenance, id="main"):
    # Make sure that the data is JSON serializable
    provenance = validate_json(provenance)

    _set_mdata_elem("provenance", value=provenance, id=id)


def add_provenance(
    loc: MuDataSlice,
    event: dict,
    id="main"
):
    provenance = get_provenance(id=id)
    provenance[loc.dehydrate()] = event
    set_provenance(provenance, id=id)


def get_mdata(full=False, id="main") -> Union[None, mu.MuData]:
    """
    Rebuild a single MuData object from all of the elements
    in the different parts of the session state.
    """

    if not get_mdata_exists(id=id):
        return

    mdata = _get_mdata_elem("data", id=id)
    assert isinstance(mdata, mu.MuData)

    if full:
        mdata.uns["mudata-explorer-process"] = get_process(id=id)
        mdata.uns["mudata-explorer-views"] = get_views(id=id)
        mdata.uns["mudata-explorer-history"] = get_history(id=id)
        mdata.uns["mudata-explorer-provenance"] = get_provenance(id=id)

    return mdata


def set_mdata(
    mdata: mu.MuData,
    timestamp: Union[None, str] = None,
    process: Union[None, str] = None,
    params: Union[None, Dict[str, Any]] = None,
    updated_keys: Union[List[str], str] = None,
    full=False,
    id="main"
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

    # Make sure that the object is of the right type
    assert isinstance(mdata, mu.MuData), type(mdata)

    # Mark that the mdata exists
    set_mdata_exists(True, id=id)

    # Save the MuData object
    _set_mdata_elem("data", value=mdata, id=id)

    if full:
        # If any of the components exist, add them to the session state
        if mdata.uns.get("mudata-explorer-views") is not None:
            set_views(mdata.uns.get("mudata-explorer-views"), id=id)
        if mdata.uns.get("mudata-explorer-process") is not None:
            set_process(mdata.uns.get("mudata-explorer-process"), id=id)
        if mdata.uns.get("mudata-explorer-settings") is not None:
            set_settings(mdata.uns.get("mudata-explorer-settings"), id=id)
        if mdata.uns.get("mudata-explorer-history") is not None:
            set_history(mdata.uns.get("mudata-explorer-history"), id=id)
        if mdata.uns.get("mudata-explorer-provenance") is not None:
            set_provenance(mdata.uns.get("mudata-explorer-provenance"), id=id)

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
        add_history(event, id=id)

        # If any updated keys were provided
        if isinstance(updated_keys, str):
            updated_keys = [updated_keys]
        if isinstance(updated_keys, list):
            for kw in updated_keys:
                add_provenance(kw, event, id=id)


def setup_mdata(id="main"):
    mdata = _setup_mdata()
    set_mdata(mdata, id=id)


def _setup_mdata() -> mu.MuData:
    mdata = mu.MuData({
        '_blank': ad.AnnData(
            X=np.array([[]]),
            obs=[],
            var=[]
        )
    })
    add_mdata_uns(mdata)
    return mdata


def list_modalities(id="main"):
    mdata: mu.MuData = get_mdata(id=id, full=False)
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


def tree_tables(axis: int, id="main") -> List[str]:
    """
    Return a list of all tables in the MuData object.
    axis = 0 for observations, 1 for variables
    """
    _is_axis(axis)

    tables = []

    if get_mdata_exists(id=id):

        # If the orientation is to observations,
        # and there is observation metadata
        if axis == 0:
            if get_mdata(id=id, full=False).obs.shape[1] > 0:
                tables.append("Observation Metadata")

        # Add tables for each modality
        tables.extend([
            join_kws(modality, table)
            for modality in list_modalities(id=id)
            for table in list_tables(modality, axis, id=id)
        ])

    return tables


def list_tables(modality: str, axis: int, id="main"):
    if not get_mdata_exists(id=id):
        return []
    _is_axis(axis)
    mdata = get_mdata(id=id, full=False)
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


def list_cnames(tables: Union[str, List[str]], axis=0, id="main"):

    _is_axis(axis)

    if not get_mdata_exists(id=id):
        return []

    mdata = get_mdata(id=id, full=False)

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
    mdata: Optional[mu.MuData] = None,
    id="main"
) -> pd.DataFrame:

    _is_axis(axis)

    # Get the complete set of data
    if mdata is None:
        mdata = get_mdata(id=id, full=False)

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
    tables: Union[List[str], str],
    axis: int,
    mdata: Optional[mu.MuData] = None,
    id="main"
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
            df = get_dataframe_table(None, "metadata", axis, mdata=mdata, id=id)
        elif '.' in table:
            modality, attr = table.split(".", 1)
            df = get_dataframe_table(modality, attr, axis, mdata=mdata, id=id)
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
    cname: str,
    id="main"
):
    if table is None:
        return

    _is_axis(axis)

    df = join_dataframe_tables(
        table,
        axis,
        mdata=mdata,
        id=id
    )
    # Make sure that the column is present
    assert cname in df.columns, f"Could not find column: {cname} in table {table} (axis={axis})"
    # Just get the column of interest
    df = df[cname]

    # If there are multiple columns with the same name
    if len(df.shape) > 1:
        df = df.iloc[:, 0]

    return df


def save_annot(
    mdata: mu.MuData,
    loc: MuDataSlice,
    column_dat: Union[pd.Series, pd.DataFrame],
    params: dict,
    process_type: str,
    figures: Optional[List[dict]],
    id="main"
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
    set_mdata(mdata, full=False, id=id)

    # Add it to the history
    add_history(event, id=id)

    # Mark the source of the table which was added
    add_provenance(loc, event, id=id)


def _overlapping_obs(mdata: mu.MuData, df: pd.DataFrame):

    # Get the observations which have been made previously
    obs = (
        set(mdata.obs.index)
        if mdata is not None
        else set()
    )

    # Add the observations from the new data
    obs.update(set(df.index))
    obs = list(obs)

    return obs


def add_modality(
    mdata: Optional[mu.MuData],
    mod_name: str,
    df: pd.DataFrame
) -> mu.MuData:
    """Add a new modality to the MuData object."""

    # If the MuData object is None, create a new one
    if mdata is None:
        mdata = setup_mdata()

    # Make sure that the modality doesn't already exist
    if isinstance(mdata, mu.MuData) and mod_name in mdata.mod:
        raise ValueError(f"Modality '{mod_name}' already exists.")

    # Add the modality name to the variable names
    # to ensure uniqueness
    df = df.rename(columns=lambda cname: f"{mod_name}:{cname}")

    # Add the new modality
    mdata.mod[mod_name] = ad.AnnData(X=df)

    # Update the total set of observation names
    mdata.update()
    return mdata


def add_obs(mdata: mu.MuData, df: pd.DataFrame):
    """
    Add a new set of observation metadata.
    This requires a dedicated method because it modifies the index
    of each of the AnnData objects.
    """

    # Get the overlap of observations between
    # the new data and the existing data
    obs = _overlapping_obs(mdata, df)

    # Update the index of the observation metadata
    df = df.reindex(index=obs)

    # Create a new MuData object
    mdata = mu.MuData(
        {
            **{
                kw: adata[obs]
                for kw, adata in (
                    mdata.mod
                    if isinstance(mdata, mu.MuData)
                    else {}
                ).items()
            }
        },
        uns=mdata.uns if mdata is not None else {}
    )
    mdata.update()

    add_mdata_uns(mdata)

    return mdata


def add_mdata_uns(mdata: mu.MuData):
    for kw, val in [
        ("views", []),
        ("process", {}),
        ("settings", {}),
        ("history", []),
        ("provenance", {})
    ]:
        if f"mudata-explorer-{kw}" not in mdata.uns:
            mdata.uns[f"mudata-explorer-{kw}"] = val


def get_supp_figs(id="main") -> List[str]:
    """Return the list of figures which are stored in the provenance."""

    return [
        f"{loc}:{ix}"
        for loc, prov in get_provenance(id=id).items()
        if prov.get("figures") is not None
        for ix, fig in enumerate(prov["figures"])
        if fig is not None
    ]
