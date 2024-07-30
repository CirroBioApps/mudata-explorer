from typing import List
from mudata_explorer import views
from mudata_explorer import process
from mudata_explorer.base.base import MuDataAppHelpers
from mudata_explorer.base.view import View
from mudata_explorer.base.process import Process
from mudata_explorer.helpers.read_table import read_table # noqa
from mudata_explorer.helpers.sanitize_types import sanitize_types # noqa
from mudata_explorer.helpers.join_kws import join_kws # noqa
import pandas as pd


def list_resources(module, iter_depth=0):
    if iter_depth > 2:
        return []
    resources = []
    for submodule in dir(module):
        if submodule.startswith("__"):
            continue
        if hasattr(getattr(module, submodule), "type"):
            if isinstance(getattr(getattr(module, submodule), "type"), str):
                resources.append(getattr(module, submodule))
        else:
            resources.extend(
                list_resources(
                    getattr(module, submodule),
                    iter_depth=iter_depth + 1
                )
            )
    return resources


all_views: List[View] = list_resources(views)
all_processes: List[Process] = list_resources(process)


def get_view_by_type(view_type: str) -> View:
    for view in all_views:
        if view.type == view_type:
            return view
    raise ValueError(f"View type '{view_type}' not found.")


def all_view_types() -> List[str]:
    return [
        view.type
        for view in all_views
        if getattr(view, "type", None) is not None
    ]


def all_process_types() -> List[str]:
    return [
        process.type
        for process in all_processes
        if getattr(process, "type", None) is not None
    ]


def make_view(type: str, **kwargs) -> View:
    view = get_view_by_type(type)
    return view(**kwargs)


def get_process_by_type(process_type: str) -> Process:
    for proc in all_processes:
        if proc.type == process_type:
            return proc
    raise ValueError(f"Process type '{process_type}' not found.")


def make_process(type: str, **kwargs) -> Process:
    proc = get_process_by_type(type)
    return proc(**kwargs)


def asset_categories(list_of_assets):

    all_categories = list(set([
        asset.category
        for asset in list_of_assets
        if hasattr(asset, "category") and getattr(asset, "category")
    ]))
    all_categories.sort()
    return all_categories


def filter_by_category(
    list_of_assets: List[MuDataAppHelpers],
    category: str
):
    return [
        asset
        for asset in list_of_assets
        if hasattr(asset, "category")
        if category == asset.category
    ]


def asset_dataframe(
    list_of_assets: List[MuDataAppHelpers]
) -> pd.DataFrame:
    """Return a DataFrame of asset types and names."""
    return pd.DataFrame(
        dict(
            name=asset.name,
            type=asset.type,
            help_text=getattr(asset, "help_text", None)
        )
        for asset in list_of_assets
    ).drop_duplicates()
