from typing import Optional
from mudata_explorer.app.mdata import get_mdata, setup_mdata
from mudata_explorer.helpers.assets import get_view_by_type
from mudata_explorer.helpers.views import get_views, set_views
import muon as mu



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


