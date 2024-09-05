from typing import Optional
from mudata_explorer.app.mdata import get_views, set_views
from mudata_explorer.helpers.assets import get_view_by_type
import muon as mu



def add_view(
    view_type: str,
    mdata: Optional[mu.MuData] = None,
    params: Optional[dict] = None,
    ix: int = -1,
    id="main"
):
    if mdata is None:
        views = get_views(id=id)
    else:
        views = mdata.uns.get("mudata-explorer-views", [])

    template = get_view_by_type(view_type).template()
    if ix == -1:
        views.append(template)
    else:
        views.insert(ix, template)
    if params is not None:
        views[ix]["params"] = params
    if mdata is None:
        set_views(views, id=id)
    else:
        mdata.uns["mudata-explorer-views"] = views
