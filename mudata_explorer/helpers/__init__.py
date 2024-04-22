from typing import List
from mudata_explorer import views
from mudata_explorer.base import View

all_views: List[View] = [
    getattr(getattr(views, view_folder), view)
    for view_folder in dir(views)
    if not view_folder.startswith("__")
    for view in dir(getattr(views, view_folder))
    if hasattr(getattr(getattr(views, view_folder), view), "mdata")
    if hasattr(getattr(getattr(views, view_folder), view), "type")
]
