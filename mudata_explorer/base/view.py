from mudata_explorer.base.base import MuDataAppHelpers
from streamlit.delta_generator import DeltaGenerator


class View(MuDataAppHelpers):

    type: str
    name: str
    help_text: str
    editable: bool
    defaults: dict = {}
    view_container: DeltaGenerator

    def __init__(
        self,
        ix: int,
        params: dict,
        editable: bool,
        **kwargs
    ):
        self.ix = ix
        self.type = type
        self.params = {
            kw: params.get(kw, val)
            for kw, val in self.get_schema_defaults(self.schema)
        }

        # Whether the params are editable
        self.params_editable = editable

    @classmethod
    def build(
        cls,
        ix=-1,
        params=dict(),
        editable=False
    ):
        return cls(
            ix=ix,
            params=params,
            editable=editable
        )

    def attach(self, container: DeltaGenerator):

        # Set up the parameters for viewing
        self.get_data(container)

        # Let the user run the method, catching any errors
        if not self.params_complete:
            container.write("Please complete all input fields")
            return

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
            help_text=cls.help_text,
            params=cls.params
        )

    def display(self):
        """Primary method which is executed to render the view."""
        pass

    def runtime_options(self, container: DeltaGenerator):
        """
        Optional method which can collection options for the view.
        The container for the view will be the sidebar which is displayed
        when the "Edit Figures" option is enabled.
        """
        pass

    def param_key(self, kw):
        return f"view-{self.ix}-{kw}"
