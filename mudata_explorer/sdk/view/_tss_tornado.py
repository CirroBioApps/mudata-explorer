# THIS FILE IS AUTOGENERATED

from mudata_explorer.sdk.helpers import collapse_params
from mudata_explorer.helpers.add_view import add_view
from mudata import MuData


def tss_tornado(mdata: MuData, **extra_params):
    """
    Visualize chromatin accessibility surrounding transcription start sites.
    """

    assert isinstance(mdata, MuData), "mdata must be a MuData object"
    extra_params = collapse_params(extra_params)

    add_view(
        'tss-tornado',
        mdata,
        params={
            
        }
    )
