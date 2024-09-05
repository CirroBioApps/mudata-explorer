from mudata_explorer.base.form import MuDataAppForm
from typing import Optional
import muon as mu


class MuDataAppAction:
    """
    Class with functions and attributes shared by Views and Processes
    """

    # Defines the position of the action within the larger app
    ix: int
    # All inputs to the action are defined as a schema
    schema: dict
    # All inputs to the action are instantiated as a form
    form: MuDataAppForm
    # Used to store unstructured data
    uns: dict
    # Unique string identifying the action
    type: str
    # Human readable name
    name: str
    # Broader category
    category: str
    # Flag used to indicate whether the parameters are editable
    params_editable: bool
    # Optional help text
    help_text: Optional[str] = False
    # Optionally attach a MuData object to the object
    mdata = Optional[mu.MuData]

    def __init__(self, ix: int):
        """Instantiate the schema as a form element."""
        self.form = MuDataAppForm(self.schema, ix=ix)
        self.ix = ix

    def display(self):
        """Overridden by child classes."""
        pass
