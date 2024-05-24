from mudata_explorer import app
from typing import List
from mudata_explorer.base.base import MuDataAppHelpers
from streamlit.delta_generator import DeltaGenerator


class View(MuDataAppHelpers):

    ix: int
    type: str
    editable: bool
    name: str
    desc: str
    categories: List[str]
    params: dict = {}
    defaults: dict = {}
    view_container: DeltaGenerator
    inputs_container: DeltaGenerator
    schema: dict = {}

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
            for kw, val in self.get_schema_defaults(self.schema)
        }

    def attach(self, container: DeltaGenerator):

        # Set up a container for the view
        self.view_container = container.container()

        if self.editable:
            # Set up an expander element for up/down/delete buttons
            self.inputs_container = container.expander("Edit Position")

        # Get the parameters from the user
        self.get_data(self.view_container)

        # Now make the display, catching any errors
        try:
            self.display(self.view_container)
        except Exception as e:
            # Log the full traceback of the exception
            self.view_container.exception(e)

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
