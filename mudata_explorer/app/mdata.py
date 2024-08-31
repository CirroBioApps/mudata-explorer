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
    mdata = _setup_mdata()
    set_mdata(mdata)


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
    mdata: mu.MuData,
    mod_name: str,
    df: pd.DataFrame
):
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


def add_history(event: dict):
    if not has_mdata():
        return
    history = get_history()
    history.insert(0, event)
    set_history(history)


def get_history(exclude=[]) -> List[dict]:
    if not has_mdata():
        return []
    else:
        return [
            h for h in json_safe(
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


def get_provenance() -> Dict[str, dict]:
    if not has_mdata():
        return {}

    mdata = get_mdata()
    assert isinstance(mdata, mu.MuData), type(mdata)

    return json_safe(mdata.uns.get("mudata-explorer-provenance", {}))


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


def get_supp_figs() -> List[str]:
    """Return the list of figures which are stored in the provenance."""

    return [
        f"{loc}:{ix}"
        for loc, prov in get_provenance().items()
        if prov.get("figures") is not None
        for ix, fig in enumerate(prov["figures"])
        if fig is not None
    ]
