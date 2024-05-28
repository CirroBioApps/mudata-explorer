from anndata import AnnData
from typing import List, Optional, Union
from muon import MuData
import pandas as pd


class MuDataSlice:
    """
    A slice of a MuData object.

    This class is used to store helper methods which work on
    a slice of a MuData object.

    The "orientation" of the slice is the axis along which the
    index is taken. If the orientation is "obs", then the index
    values of the DataFrame which represents the data will contain
    elements which are the observations.

    For `attr` and `subattr`, recognized values are:
        - None: No attribute or subattribute is used.
        - str: The name of an attribute or subattribute.
        - List[str]: A list of names of the attributes or subattributes.

    Example 1:
        modality:       None
        orientation:    "obs"
        slot:           "obs"
        attr:           None
        subattr:        None

        Resolves to: mdata.obs

    Example 2:
        modality:       None
        orientation:    "var"
        slot:           "obs"
        attr:           None
        subattr:        None

        Resolves to: mdata.obs.T

    Example 3:
        modality:       "rna"
        orientation:    "obs"
        slot:           "obsm"
        attr:           "X_pca"
        subattr:        "PC1"

        Resolves to: mdata.obsm["X_pca"]["PC1"]

    """

    def __init__(
        self,
        mdata: MuData,
        slot: str,
        orientation: str = "obs",
        modality: Optional[str] = None,
        attr: Optional[Union[str, List[str]]] = None,
        subattr: Optional[Union[str, List[str]]] = None
    ):
        self.mdata = mdata
        self.slot = slot
        self.orientation = orientation
        self.modality = modality
        self.attr = attr
        self.subattr = subattr

        assert isinstance(self.mdata, MuData)
        assert isinstance(slot, str)
        assert orientation in ["obs", "var"], f"Unexpected: {orientation}"

        if self.slot != "X":
            msg = f"Cannot access slot {slot} with orientation {orientation}"
            assert self.slot.startswith(self.orientation), msg

        if self.modality is None:
            msg = "If the modality is None, the slot must be 'obs'"
            assert self.slot == "obs", msg
        else:
            assert self.modality in self.mdata.mod.keys()

    def dataframe(self) -> pd.DataFrame:
        """Return a DataFrame containing the data in the slice."""

        if self.modality is None:
            dat: MuData = self.mdata
        else:
            dat: AnnData = self.mdata.mod[self.modality]

        if self.slot == "X":
            dat = dat.to_df()
            if self.orientation == "var":
                dat = dat.T
        else:
            assert hasattr(dat, self.slot)
            dat = getattr(dat, self.slot)

        if self.attr is not None:
            dat = dat[self.attr]

        if self.subattr is not None:
            dat = dat[self.subattr]

        return dat
