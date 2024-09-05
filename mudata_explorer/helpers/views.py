from mudata_explorer.app.mdata import get_views, set_views


def delete_view(ix: int, id="main"):
    views = get_views(id=id)
    views.pop(ix)
    set_views(views, id=id)


def duplicate_view(ix: int, id="main"):
    views = get_views(id=id)
    views.insert(ix, views[ix])
    set_views(views, id=id)
