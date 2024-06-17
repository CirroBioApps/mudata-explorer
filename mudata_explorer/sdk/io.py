from anndata import AnnData
import logging
from mudata_explorer.helpers.mudata import add_mdata_uns
from mudata_explorer.app import get_dat_hash
from muon import MuData
from pandas import DataFrame
from typing import Dict

# Set up logging
logging.basicConfig(level=logging.INFO)


def build_mdata(dat: Dict[str, DataFrame], obs=DataFrame) -> MuData:
    """Build a MuData object from a dictionary of AnnData objects."""
    mdata = MuData(
        {
            kw: AnnData(
                X=df.rename(columns=lambda cname: f"{kw}:{cname}")
            )
            for kw, df in dat.items()
        }
    )
    mdata.obs = obs
    add_mdata_uns(mdata)

    return mdata


def write_h5mu(mdata: MuData, prefix: str):
    """
    Save the MuData object to an .h5mu file using the filename:
    {prefix}-{hash}.h5mu
    """

    dat, hash, size = get_dat_hash(mdata)

    with open(f"{prefix}-{hash}.h5mu", "wb") as f:
        f.write(dat)

    logging.info(f"Saved MuData object to {prefix}-{hash}.h5mu ({size} bytes)")


def write_zarr(mdata: MuData, prefix: str):
    """
    Save the MuData object to an .zarr file using the filename:
    {prefix}-{hash}.zarr
    """

    dat, hash, size = get_dat_hash(mdata)

    dat.write_zarr(f"{prefix}-{hash}.zarr")

    logging.info(f"Saved MuData object to {prefix}-{hash}.zarr")
