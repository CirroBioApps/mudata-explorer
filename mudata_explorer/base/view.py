from typing import List
from mudata_explorer.base.base import MuDataAppHelpers
from mudata_explorer.app.mdata import set_mdata
import streamlit as st
from streamlit.delta_generator import DeltaGenerator


class View(MuDataAppHelpers):

    type: str
    name: str
    help_text: str
    editable: bool
    defaults: dict = {}

    def __init__(
        self,
        ix: int,
        params: dict,
        editable: bool,
        uns={},
        **kwargs
    ):
        self.ix = ix
        self.type = type
        self.params = {
            kw: params.get(kw, val)
            for kw, val in self.get_schema_defaults(self.schema)
        }
        self.uns = uns

        # Whether the params are editable
        self.params_editable = editable

    @classmethod
    def build(
        cls,
        ix=-1,
        params=dict(),
        uns=dict(),
        editable=False
    ):
        return cls(
            ix=ix,
            params=params,
            uns=uns,
            editable=editable
        )

    def attach(self):

        # Set up the parameters for viewing
        self.get_data()

        # Let the user run the method, catching any errors
        if not self.params_complete:
            st.write("Please complete all input fields")
            return

        # Now make the display, catching any errors
        try:
            self.display()
        except Exception as e:
            # Log the full traceback of the exception
            st.exception(e)

    @classmethod
    def template(cls):
        return dict(
            type=cls.type,
            name=cls.name,
            help_text=cls.help_text,
            params=cls.params
        )

    def display(self):
        """Primary method which is executed to render the view."""
        pass

    def runtime_options(self):
        """
        Optional method which can collection options for the view.
        The container for the view will be the sidebar which is displayed
        when the "Edit Figures" option is enabled.
        """
        pass

    def param_key(self, kw):
        return f"view-{self.ix}-{kw}"
