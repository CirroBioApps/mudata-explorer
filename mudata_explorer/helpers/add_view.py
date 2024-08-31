from typing import Optional
from mudata_explorer.app.mdata import get_mdata, setup_mdata
from mudata_explorer.helpers.assets import get_view_by_type
from mudata_explorer.helpers.views import get_views, set_views
import muon as mu



def add_view(
    view_type: str,
    mdata: Optional[mu.MuData] = None,
    params: Optional[dict] = None,
    ix: int = -1
):
    if mdata is None:
        if get_mdata() is None:
            setup_mdata()
    views = get_views(mdata)
    template = get_view_by_type(view_type).template()
    if ix == -1:
        views.append(template)
    else:
        views.insert(ix, template)
    if params is not None:
        views[ix]["params"] = params
    set_views(views, mdata)
