from typing import Dict, List, Optional, Union, Generator
import numpy as np
import pandas as pd
from mudata_explorer.app.mdata import get_mdata
from mudata_explorer.app.mdata import save_annot
from mudata_explorer.app.mdata import get_process, set_process
from mudata_explorer.base.base import MuDataAppAction
from mudata_explorer.base.slice import MuDataSlice
from mudata_explorer.helpers.params import nest_params
from streamlit.delta_generator import DeltaGenerator
from mudata import MuData


class Process(MuDataAppAction):

    ix = -1
    outputs = Dict[str, dict]

    def __init__(
        self,
        params: dict = {},
        mdata: Optional[MuData] = None,
        params_editable=True
    ):
        """
        If the optional mdata is provided, the process will
        save results to that object instead of the global object.
        """

        # Instantiate the self.form attribute from the schema
        super().__init__(-1)
        # Load any params
        self.form.load(nest_params(params))
        self.params_editable = params_editable
        self.mdata = mdata

    def get_output_locs(self) -> List[MuDataSlice]:
        """
        Return the list of output locations for this process.
        """

        # Resolve each of the output locations which have been defined
        # Each output schema may resolve to multiple locations
        # (e.g.), if multiple modalities are selected with axis=1
        locs = []
        for output in self.outputs.values():
            for loc in self.resolve_output_loc(output):
                locs.append(loc)
        return locs

    def resolve_output_loc(
        self,
        output: dict
    ) -> Generator[MuDataSlice, None, None]:

        # If the axis value ends with .T, it means that the axis
        # is transposed
        transposed = output["axis"].endswith(".T")

        # First, resolve any of the elements of the output
        # which refer to keywords in the params scope
        output = {
            kw: (
                self.params[val]
                if val in self.params
                else val
            )
            for kw, val in output.items()
        }

        if isinstance(output["modality"], np.ndarray):
            output["modality"] = output["modality"].tolist()

        # If the axis should be transposed
        if transposed:
            output["axis"] = self.params[output["axis"][:-2]]
            output["axis"] = 1 if output["axis"] == 0 else 0

        # The output type can be either a Series or DataFrame
        assert output["type"] in (pd.Series, pd.DataFrame)

        # If the output slot is None, resolve it depending on the axis
        # and the data type
        # DataFrames are put into obsm/varm
        # Series are put into obs/var
        if output.get("slot") is None:
            if output["axis"] == 0:
                if output["type"] == pd.Series:
                    # Observation metadata is stored in the .obs slot
                    output["modality"] = None
                    output["slot"] = "obs"
                else:
                    output["slot"] = "obsm"
            else:
                output["slot"] = (
                    "var" if output["type"] == pd.Series
                    else "varm"
                )

        # If the modality is an empty list,
        # it means that the output cannot be set
        if output["modality"] == []:
            return

        # If the modality is a list
        if isinstance(output["modality"], list):

            # Parse out the modality from each string
            output["modality"] = list(set([
                mod.split(".")[0]
                for mod in output["modality"]
            ]))

            # Return multiple locations for each unique modality
            for mod in output["modality"]:
                yield MuDataSlice(
                    axis=output["axis"],
                    modality=mod,
                    slot=output["slot"],
                    attr=output["attr"]
                )

        # If the modality is a single string or None
        else:
            assert isinstance(
                output["modality"],
                (list, str, type(None))
            ), type(output["modality"])

            yield MuDataSlice(
                axis=output["axis"],
                modality=(
                    output["modality"]
                    if output["modality"] is None
                    else
                    output["modality"].split(".")[0]
                ),
                slot=output["slot"],
                attr=output["attr"]
            )

    def run(self, container: DeltaGenerator):

        pass

    def save_results(
        self,
        output_kw: str,
        res: Union[pd.Series, pd.DataFrame],
        figures: Optional[List[dict]] = None,
        id="main"
    ):

        assert output_kw in self.outputs, f"Output {output_kw} not found"

        # Get the MuData object
        if self.mdata is None:
            mdata = get_mdata(id=id, full=False)
        else:
            mdata = self.mdata

        # Get the output location
        for loc in self.resolve_output_loc(self.outputs[output_kw]):
            assert isinstance(loc, MuDataSlice)

            # Save the results to the MuData object
            save_annot(
                mdata,
                loc,
                res,
                self.dehydrate(),
                self.type,
                figures,
                id=id
            )

    def execute(self):
        pass

    def param_key(self, kw):
        return f"process-{kw}"

    def update_view_param(self, kw, value, id="main"):
        # Get the current process
        process = get_process(id=id)

        # Modify the value of this param for this view
        if "params" not in process:
            process["params"] = {}

        # Set the value of the parameter
        process["params"][kw] = value

        # Save the process definition
        set_process(process=process, id=id)

        # Also update the params object
        self.params[kw] = value

    def dehydrate(self):
        """Only save those params which can be loaded."""

        return self.form.dehydrate()

    @classmethod
    def hydrate(cls, params: dict):
        return cls(params)

    def populate_params(self):
        """Populate the params attribute of the object."""

        self.params = self.form.dump(self.mdata)
