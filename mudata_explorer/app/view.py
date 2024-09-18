from mudata_explorer.app.mdata import get_views
from mudata_explorer.helpers.assets import make_view


def make_views():

    views = get_views()

    return [
        make_view(
            ix=ix,
            type=view["type"],
            params=view["params"]
        )
        for ix, view in enumerate(views)
    ]


def view_non_editable():

    # All of the views defined in the dataset
    mdata_views = make_views()

    for view in mdata_views:

        # Attach the view to the display
        view.attach()