import json
from anndata import AnnData
from typing import List, Optional, Union
from muon import MuData
import numpy as np
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

        # If the slot is .obs, the modality will be coerced to None
        if self.slot == "obs":
            self.modality = None

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
        address = (
            f"{self.modality}.{self.slot}"
            if self.modality is not None
            else self.slot
        )

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

        if mdata is None:
            return None

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

        assert isinstance(dat, (pd.Series, pd.DataFrame))

        dat.index.name = self.orientation

        return dat.dropna(how="all")

    def write(self, mdata: MuData, dat: Union[pd.Series, pd.DataFrame]):
        """Write the data in the slice to the container."""

        assert isinstance(dat, (pd.Series, pd.DataFrame)), type(dat)

        # Writing to the observation metadata
        # Note that this is a special case
        if self.slot == "obs":
            assert self.subattr is None
            assert self.orientation == "obs"

            self._write_obs(mdata, dat)

        # Writing to the variable metadata
        elif self.slot == "var":

            # The modality must be defined
            assert self.modality is not None
            assert self.subattr is None
            assert self.orientation == "var"

            # Set the index to the variables
            dat = dat.reindex(mdata.mod[self.modality].var.index)

            # If writing a column
            if self.attr is not None:
                assert isinstance(dat, pd.Series)
                mdata.mod[self.modality].var[self.attr] = dat

            # If writing a DataFrame
            else:
                mdata.mod[self.modality].var = dat

        # Writing to obsm/obsp/varm/varp
        elif self.slot in ["obsm", "obsp", "varm", "varp"]:

            # The correct orientation must be used
            assert self.slot.startswith(self.orientation)

            # If writing a DataFrame
            if isinstance(dat, pd.DataFrame):
                # No subattribute is defined
                assert self.subattr is None

                # Align the axis to the modality
                index = getattr(
                    mdata.mod[self.modality],
                    f"{self.orientation}_names"
                )
                dat = dat.reindex(index=index)
                print(dat)
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

    def _write_obs(self, mdata: MuData, dat: Union[pd.Series, pd.DataFrame]):
        """
        Writing to the .obs is a special case.

        The index of dat may contain observations which are not in the
        current mdata object.

        To ensure that all values may be added, a dummy modality will be
        added which contains the observations which are in the provided
        metadata table.
        """

        # If writing a column, simply add to the existing table
        if self.attr is not None:
            assert isinstance(dat, pd.Series)
            mdata.obs[self.attr] = dat

        # If writing a DataFrame
        else:

            # If any of the observations are not in the current mdata object
            if len(set(dat.index) - set(mdata.obs.index)) > 0:

                # Add a new modality
                mdata.mod["_obs"] = AnnData(
                    X=pd.DataFrame(
                        np.zeros((dat.shape[0], 1)),
                        index=dat.index
                    )
                )

                # Update the complete set of observations
                mdata.update()

            # Set the metadata
            mdata.obs = dat.reindex(
                index=mdata.obs_names
            )
