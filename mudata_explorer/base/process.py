from typing import List, Optional, Tuple, Union
import pandas as pd
import muon as mu
from mudata_explorer import app
from mudata_explorer.base.base import MuDataAppHelpers
from mudata_explorer.base.slice import MuDataSlice
from scipy.stats import zscore
from streamlit.delta_generator import DeltaGenerator


class Process(MuDataAppHelpers):

    type: str
    name: str
    desc: str
    categories: List[str]
    schema: dict
    ix = -1
    output_type: Union[pd.Series, pd.DataFrame]
    output_modalities: List[str]
    figures: Optional[List[dict]] = None

    def __init__(
        self,
        params: dict = {}
    ):
        self.params = {
            kw: params.get(kw, val)
            for kw, val in self.get_schema_defaults(self.schema)
        }

    def run(self, container: DeltaGenerator):

        pass

    def get_data(self, container: DeltaGenerator):
        """The difference between .get_data for the Process and View
        is that the Process needs to know where to store the results."""

        # Run the method of the parent class
        super().get_data(container)

        # Now figure out what the output modalities are

        # If the orientation is to obs
        if self.orientation == "observations":

            # Ask the user which modality to save the outputs to
            all_modalities = app.list_modalities()
            dest_modality = container.selectbox(
                "Write results to modality:",
                all_modalities
            )
            self.output_modalities = [dest_modality]

        # If the orientation is to var
        else:

            # Look through the tables which were selected by the user
            self.output_modalities = self._find_modalities()

        # Make sure that all of the output modalities are valid
        assert all([
            mod in app.list_modalities()
            for mod in self.output_modalities
        ]), self.output_modalities

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
