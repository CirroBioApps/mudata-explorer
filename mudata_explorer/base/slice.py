import json
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

    orientation: str
    modality: str
    slot: str
    attr: Optional[str]
    subattr: Optional[str]
    _allowed_slots = ["obs", "var", "X", "obsm", "obsp", "varm", "varp"]

    def __init__(
        self,
        slot: str,
        orientation: str = "obs",
        modality: Optional[str] = None,
        attr: Optional[Union[str, List[str]]] = None,
        subattr: Optional[Union[str, List[str]]] = None
    ):
        assert slot in self._allowed_slots, slot
        self.slot = slot
        self.orientation = orientation
        self.modality = modality
        self.attr = attr
        self.subattr = subattr

        assert isinstance(slot, str)
        assert orientation in ["obs", "var"], f"Unexpected: {orientation}"

        if self.slot != "X":
            msg = f"Cannot access slot {slot} with orientation {orientation}"
            assert self.slot.startswith(self.orientation), msg

        if self.modality is None:
            msg = "If the modality is None, the slot must be 'obs'"
            assert self.slot == "obs", msg

    @property
    def params(self):
        return dict(
            slot=self.slot,
            orientation=self.orientation,
            modality=self.modality,
            attr=self.attr,
            subattr=self.subattr
        )

    def dehydrate(self):
        return json.dumps(self.params, sort_keys=True)

    @classmethod
    def hydrate(cls, params: str):
        return cls(**json.loads(params))

    @property
    def address(self) -> str:
        """Return a string representation of the slice."""
        address = f"{self.modality}.{self.slot}"

        for key in ["attr", "subattr"]:
            if getattr(self, key) is not None:
                vals = getattr(self, key)
                if isinstance(vals, list):
                    address += f"[{', '.join(vals)}]"
                else:
                    address += f"[{vals}]"

        return address

    def dataframe(self, mdata: MuData) -> pd.DataFrame:
        """Return a DataFrame containing the data in the slice."""

        assert isinstance(mdata, MuData)

        if self.modality is None:
            dat: MuData = mdata
        else:
            dat: AnnData = mdata.mod[self.modality]

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

    def write(self, mdata: MuData, dat: Union[pd.Series, pd.DataFrame]):
        """Write the data in the slice to the container."""

        assert self.attr is not None

        # Writing to the observation metadata
        if self.slot == "obs":
            # Can only write a column
            assert isinstance(dat, pd.Series)
            assert self.subattr is None
            assert self.orientation == "obs"
            mdata.obs[self.attr] = dat

        # Writing to the variable metadata
        elif self.slot == "var":
            # Can only write a column
            assert isinstance(dat, pd.Series)
            assert self.subattr is None
            assert self.orientation == "var"
            mdata.var[self.attr] = dat

        # Writing to obsm/obsp/varm/varp
        elif self.slot in ["obsm", "obsp", "varm", "varp"]:

            # The correct orientation must be used
            assert self.slot.startswith(self.orientation)

            # If writing a DataFrame
            if isinstance(dat, pd.DataFrame):
                # No subattribute is defined
                assert self.subattr is None
                getattr(mdata.mod[self.modality], self.slot)[self.attr] = dat
            # If writing a Series
            else:
                # A subattribute is defined
                assert self.subattr is not None
                (
                    getattr(
                        mdata.mod[self.modality],
                        self.slot
                    )
                    [self.attr]
                    [self.subattr]
                ) = dat
