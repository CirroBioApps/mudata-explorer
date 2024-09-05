import hashlib
from typing import Union, Optional
from mudata_explorer.app.mdata import add_history
from mudata_explorer.helpers.timestamp import get_timestamp
from mudata_explorer.app.mdata import get_mdata, get_mdata_exists
import muon as mu
from mudata_explorer.helpers.io import mdata_to_binary


def hash_dat(dat, n: Union[int, None] = 16):
    """Compute the hash of an object."""
    hash = hashlib.sha256()
    hash.update(dat)
    hex = hash.hexdigest()
    if n is not None:
        hex = hex[:n]
    return hex


def set_mdata_hash(mdata_hash: str, id="main"):
    """Record the hash of the data in the history."""
    add_history(dict(
        process="data_hash",
        params=dict(
            hash=mdata_hash
        ),
        timestamp=get_timestamp(),
        id=id
    ))


def get_dat_hash(mdata: Optional[mu.MuData] = None, id="main"):
    if mdata is None and get_mdata_exists(id=id) is False:
        return None, None, None

    if mdata is None:
        mdata = get_mdata(full=True, id=id)

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
