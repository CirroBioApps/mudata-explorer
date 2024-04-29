from typing import List
from mudata_explorer import views
from mudata_explorer.base.view import View
from mudata_explorer.helpers.read_table import read_table # noqa
from mudata_explorer.helpers.sanitize_types import sanitize_types # noqa

all_views: List[View] = [
    getattr(getattr(views, view_folder), view)
    for view_folder in dir(views)
    if not view_folder.startswith("__")
    for view in dir(getattr(views, view_folder))
    if hasattr(getattr(getattr(views, view_folder), view), "type")
]


def get_view_by_type(view_type: str) -> View:
    for view in all_views:
        if view.type == view_type:
            return view
    raise ValueError(f"View type '{view_type}' not found.")


def make_view(type: str, **kwargs) -> View:
    view = get_view_by_type(type)
    return view(type=type, **kwargs)
