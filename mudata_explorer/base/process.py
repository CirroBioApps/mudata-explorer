from typing import List, Optional, Union
import pandas as pd
from mudata_explorer import app
from mudata_explorer.base.base import MuDataAppHelpers
from mudata_explorer.base.slice import MuDataSlice
from streamlit.delta_generator import DeltaGenerator


class Process(MuDataAppHelpers):

    ix = -1
    output_type: Union[pd.Series, pd.DataFrame]
    figures: Optional[List[dict]] = None

    def __init__(
        self,
        params: dict = {}
    ):
        self.params = {
            kw: params.get(kw, val)
            for kw, val in self.get_schema_defaults(self.schema)
        }
        # Params are always editable for a new process
        self.params_editable = True

    def get_output_locs(self, dest_key) -> List[MuDataSlice]:
        """
        Return the list of output locations for this process.
        """

        # Look through the tables which were selected by the user
        # in order to find the modalities used for inputs
        output_modalities = self._find_modalities()

        if len(output_modalities) == 0:
            return []

        # If the orientation is to the observations
        if self.orientation == "observations":

            # If the output is a Series (a single column)
            if self.output_type == pd.Series:

                # Write results to the mdata.obs slot
                return [
                    MuDataSlice(
                        orientation="obs",
                        modality=None,
                        slot="obs",
                        attr=dest_key
                    )
                ]

            # If the output is a DataFrame
            else:

                # Write results to the .obsm slot of the
                # first modality used in the process
                return [
                    MuDataSlice(
                        orientation="obs",
                        modality=output_modalities[0],
                        slot="obsm",
                        attr=dest_key
                    )
                ]
        # If the orientation is to variables
        else:

            # If the output is a Series (a single row)
            if self.output_type == pd.Series:
                # Save to .var
                slot = "var"

            # If the output is a DataFrame
            else:
                # Save to .varm
                slot = "varm"

            # Save results in all of the modalities used in the process
            return [
                MuDataSlice(
                    orientation="var",
                    modality=modality,
                    slot=slot,
                    attr=dest_key
                )
                for modality in output_modalities
            ]

    def run(self, container: DeltaGenerator):

        pass

    def _find_modalities(self) -> List[str]:
        """Look through all of the params and return a list of all
        of the modalities which have been selected for source data."""

        return list(set(
            [
                # If the user selects a table
                table.split(".", 1)[0]
                for kw, vals in self.params.items()
                if kw.endswith('.tables')
                for table in vals
            ] + [
                # If the user selects a column
                val.split(".", 1)[0]
                for kw, val in self.params.items()
                if kw.endswith('.modality')
            ]
        ))

    def execute(self) -> Union[pd.Series, pd.DataFrame]:
        pass

    def param_key(self, kw):
        return f"process-{kw}"

    def update_view_param(self, kw, value):
        # Get the MuData object
        mdata = app.get_mdata()

        # Modify the value of this param for this view
        if "params" not in mdata.uns["mudata-explorer-process"]:
            mdata.uns["mudata-explorer-process"]["params"] = {}

        # Set the value of the parameter
        mdata.uns["mudata-explorer-process"]["params"][kw] = value

        # Save the MuData object
        app.set_mdata(mdata)

        # Also update the params object
        self.params[kw] = value

    def locate_results(
        self,
        dest_modality: str,
        dest_key: str
    ) -> MuDataSlice:
        """Determine the location where a set of results will be saved."""

        # Depending on the orientation, set the destination attribute
        if self.params["orientation"] == "observations":
            attr = "obs"
        else:
            attr = "var"

        # DataFrames get written as their own table
        if self.output_type == pd.DataFrame:
            attr = attr + "m"

        # Return the location which was written
        return MuDataSlice(
            orientation=self.params["orientation"][:3],
            modality=dest_modality,
            slot=attr,
            attr=dest_key
        )

    def save_results(
        self,
        loc: MuDataSlice,
        res: Union[pd.Series, pd.DataFrame]
    ) -> MuDataSlice:

        # Get the MuData object
        mdata = app.get_mdata()

        # Save the results to the MuData object
        app.save_annot(
            mdata,
            loc,
            res,
            self.dehydrate(),
            self.type,
            self.figures
        )

    def dehydrate(self):
        """Only save those params which can be loaded."""

        return {
            kw: self.params[kw]
            for kw, _ in self.get_schema_defaults(self.schema)
        }

    @classmethod
    def hydrate(cls, params: dict):
        return cls(params)
