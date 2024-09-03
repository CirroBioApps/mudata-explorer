from mudata_explorer.app.mdata import get_mdata, set_mdata, has_mdata
from mudata_explorer.app.mdata import get_supp_figs
from mudata_explorer.app.mdata import tree_tables, list_cnames, join_dataframe_tables, get_dataframe_column
from mudata_explorer.app.query_params import get_editable_flag
from mudata_explorer.base import all_transforms, get_transform
from mudata_explorer.base.form import MuDataAppForm
from mudata_explorer.helpers.join_kws import join_kws
from typing import Any, Dict, List, Optional
import muon as mu
import pandas as pd
import plotly.express as px
import streamlit as st
from streamlit.delta_generator import DeltaGenerator


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
