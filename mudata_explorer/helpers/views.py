from typing import Optional, List
from mudata_explorer.helpers.io import validate_json, json_safe
from mudata_explorer.app.mdata import get_mdata, set_mdata
import muon as mu


def delete_view(ix: int):
    views = get_views()
    views.pop(ix)
    set_views(views)


def duplicate_view(ix: int):
    views = get_views()
    views.insert(ix, views[ix])
    set_views(views)


def get_views(mdata: Optional[mu.MuData] = None) -> List[dict]:
    if mdata is None:
        mdata = get_mdata()
    if mdata is None:
        return []
    assert isinstance(mdata, mu.MuData), type(mdata)
    return json_safe(mdata.uns.get("mudata-explorer-views", []))


def set_views(views, mdata: Optional[mu.MuData] = None):
    use_global = mdata is None
    if use_global:
        mdata = get_mdata()

    # Make sure that the data is JSON serializable
    views = validate_json(views)

    mdata.uns["mudata-explorer-views"] = views
    if use_global:
        set_mdata(mdata)