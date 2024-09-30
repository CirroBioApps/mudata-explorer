from copy import copy
import json
from pathlib import Path
import muon as mu
from mudata_explorer.app.mdata import set_mdata
from functools import lru_cache


@lru_cache
def _read_local_file(filename: str) -> mu.MuData:
    file = Path(__file__).parent / filename
    mdata = mu.read_h5mu(file)
    # Parse the mudata elements as JSON strings
    for suffix in ["history", "views", "provenence", "settings", "history"]:
        kw = f"mudata-explorer-{suffix}"
        if kw in mdata.uns:
            if isinstance(mdata.uns[kw], str):
                mdata.uns[kw] = json.loads(mdata.uns[kw])
    return mdata


def load_explanation(filename="microbiome-report-176ff7077db10524.h5mu"):
    """
    Read in the microbiome-report.h5mu file in this directory
    and save it to the session state.
    """
    mdata = _read_local_file(filename)
    set_mdata(copy(mdata), full=True)
