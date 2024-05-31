from anndata import AnnData
from muon import MuData
import numpy as np
import pandas as pd


def _overlapping_obs(mdata: MuData, df: pd.DataFrame):

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
    mdata: MuData,
    mod_name: str,
    df: pd.DataFrame
):
    """Add a new modality to the MuData object."""

    # If the MuData object is None, create a new one
    if mdata is None:
        mdata = setup_mdata()

    # Make sure that the modality doesn't already exist
    if isinstance(mdata, MuData) and mod_name in mdata.mod:
        raise ValueError(f"Modality '{mod_name}' already exists.")

    # Add the modality name to the variable names
    # to ensure uniqueness
    df = df.rename(columns=lambda cname: f"{mod_name}:{cname}")

    # Add the new modality
    mdata.mod[mod_name] = AnnData(X=df)

    # Update the total set of observation names
    mdata.update()
    return mdata


def add_obs(mdata: MuData, df: pd.DataFrame):
    """
    Add a new set of observation metadata.
    This requires a dedicated method because it modifies the index
    of each of the AnnData objects.
    """

    # Get the overlap of observations between the new data and the existing data
    obs = _overlapping_obs(mdata, df)

    # Update the index of the observation metadata
    df = df.reindex(index=obs)

    # Create a new MuData object
    mdata = MuData(
        {
            **{
                kw: adata[obs]
                for kw, adata in (
                    mdata.mod
                    if isinstance(mdata, MuData)
                    else {}
                ).items()
            }
        },
        uns=mdata.uns if mdata is not None else {}
    )
    mdata.update()

    add_mdata_uns(mdata)

    return mdata


def setup_mdata() -> MuData:
    mdata = MuData({
        '_blank': AnnData(
            X=np.array([[]]),
            obs=[],
            var=[]
        )
    })
    add_mdata_uns(mdata)
    return mdata


def add_mdata_uns(mdata: MuData):
    for kw, val in [
        ("views", []),
        ("process", {}),
        ("settings", {}),
        ("history", []),
        ("provenance", {})
    ]:
        if f"mudata-explorer-{kw}" not in mdata.uns:
            mdata.uns[f"mudata-explorer-{kw}"] = val
