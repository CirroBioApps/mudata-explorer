import streamlit as st
from typing import List
from mudata_explorer.base.base import MuDataAppHelpers
from streamlit.delta_generator import DeltaGenerator


class View(MuDataAppHelpers):

    ix: int
    type: str
    editable: bool
    params: dict
    name: str
    desc: str
    categories: List[str]
    params: dict = {}
    defaults: dict = {}
    view_container: DeltaGenerator
    inputs_container: DeltaGenerator

    def __init__(
        self,
        ix: int,
        type: str,
        editable: bool,
        name: str,
        desc: str,
        params: dict
    ):
        self.ix = ix
        self.type = type
        self.editable = editable
        self.name = name
        self.desc = desc
        self.params = {
            kw: params.get(kw, val)
            for kw, val in self.defaults.items()
        }

    def attach(self, container: DeltaGenerator):

        # Set up a container for the view
        self.view_container = container.container()

        if self.editable:
            # Set up an expander element for the parameters
            self.inputs_container = container.expander("Edit Settings")

        # Now make the display
        self.display(self.view_container)

    @classmethod
    def template(cls):
        return dict(
            type=cls.type,
            name=cls.name,
            desc=cls.desc,
            params=cls.params
        )

    def display(self):
        """Primary method which is executed to render the view."""
        pass

    def param_key(self, kw):
        return f"view-{self.ix}-{kw}"

    def input_value_kwargs(self, kw):
        """Each input value element will be populated with default kwargs."""
        return dict(
            value=self.params[kw],
            key=self.param_key(kw),
            on_change=self.input_value_change,
            args=(kw,)
        )
    
    def input_value_change(self, kw):
        # Get the value provided by the user
        value = st.session_state[self.param_key(kw)]

        # If this is different than the params
        if value != self.params[kw]:

            # Update the view in the mdata object
            self.update_view_param(kw, value)

    def update_view_param(self, kw, value):
        # Get the MuData object
        mdata = self.get_mdata()

        # Modify the value of this param for this view
        mdata.uns["mudata-explorer-views"][self.ix]["params"][kw] = value

        # Save the MuData object
        self.set_mdata(mdata)
