from typing import List
from mudata_explorer import views
from mudata_explorer import process
from mudata_explorer.base.view import View
from mudata_explorer.base.process import Process
from mudata_explorer.helpers.read_table import read_table # noqa
from mudata_explorer.helpers.sanitize_types import sanitize_types # noqa


def list_resources(module):
    return [
        getattr(getattr(module, view_folder), view)
        for view_folder in dir(module)
        if not view_folder.startswith("__")
        for view in dir(getattr(module, view_folder))
        if hasattr(getattr(getattr(module, view_folder), view), "type")
    ]


all_views: List[View] = list_resources(views)
all_processes: List[Process] = list_resources(process)


def get_view_by_type(view_type: str) -> View:
    for view in all_views:
        if view.type == view_type:
            return view
    raise ValueError(f"View type '{view_type}' not found.")


def make_view(type: str, editable: bool, **kwargs) -> View:
    view = get_view_by_type(type)
    return view(type=type, editable=editable, **kwargs)


def get_process_by_type(process_type: str) -> Process:
    for proc in all_processes:
        if proc.type == process_type:
            return proc
    raise ValueError(f"Process type '{process_type}' not found.")


def make_process(type: str, **kwargs) -> Process:
    proc = get_process_by_type(type)
    return proc(type=type, **kwargs)


def asset_categories(list_of_assets):

    all_categories = list(set([
        category
        for asset in list_of_assets
        for category in asset.categories
    ]))
    all_categories.sort()
    return all_categories


def filter_by_category(list_of_assets, category):
    return [
        asset
        for asset in list_of_assets
        if category in asset.categories
    ]


def asset_name_list(list_of_assets):
    return [view.name for view in list_of_assets]


def asset_type_list(list_of_assets):
    return [view.type for view in list_of_assets]


def asset_desc_list(list_of_assets):
    return [f"{asset.name}: {asset.desc}" for asset in list_of_assets]


def asset_type_desc_lists(list_of_assets):
    return (
        asset_type_list(list_of_assets),
        asset_desc_list(list_of_assets)
    )
