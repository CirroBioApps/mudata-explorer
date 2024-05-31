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

        # Whether the params are editable is driven by the settings
        settings = app.get_settings() if app.has_mdata() else {}
        self.params_editable = settings.get("editable", True)

    def attach(self, container: DeltaGenerator):

        # Set up a container for the view
        self.view_container = container.container()

        if self.editable:
            # Set up an expander element for up/down/delete buttons
            self.inputs_container = container.expander("Edit Display")

        # Get the parameters from the user
        self.get_data(self.inputs_container)

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
