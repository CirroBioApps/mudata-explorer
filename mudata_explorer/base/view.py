from mudata_explorer.base.base import MuDataAppHelpers
from streamlit.delta_generator import DeltaGenerator


class View(MuDataAppHelpers):

    editable: bool
    defaults: dict = {}
    view_container: DeltaGenerator

    def __init__(
        self,
        ix: int,
        type: str,
        name: str,
        desc: str,
        params: dict,
        editable: bool
    ):
        self.ix = ix
        self.type = type
        self.name = name
        self.desc = desc
        self.params = {
            kw: params.get(kw, val)
            for kw, val in self.get_schema_defaults(self.schema)
        }

        # Whether the params are editable
        self.params_editable = editable

    def attach(self, container: DeltaGenerator):

        # Set up the parameters for viewing
        self.get_data(container)

        # Now make the display, catching any errors
        try:
            self.display(container)
        except Exception as e:
            # Log the full traceback of the exception
            container.exception(e)

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
