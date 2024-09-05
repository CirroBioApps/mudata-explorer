from typing import Optional
from mudata_explorer.app.mdata import get_view, set_view
from mudata_explorer.helpers.params import nest_params
from mudata_explorer.base.base import MuDataAppAction
import streamlit as st
from streamlit.delta_generator import DeltaGenerator


class View(MuDataAppAction):

    type: str
    name: str
    help_text: str
    editable: bool
    defaults: dict = {}
    params = {}

    def __init__(
        self,
        ix: int,
        params: dict,
        uns={}
    ):
        # Instantiate the self.form attribute from the schema
        # (also set self.ix)
        super().__init__(ix)

        # Load any params
        self.form.load(nest_params(params))
        self.params = self.form.dehydrate()

        # Set any unstructured data
        self.uns = uns

    @classmethod
    def build(
        cls,
        ix=-1,
        params=dict(),
        uns=dict()
    ):
        return cls(
            ix=ix,
            params=params,
            uns=uns
        )
    
    def display_form(self):
        """Render the form element, and then save changes."""
        self.form.render()
        self.save_changes()
        # Save the output of the form as self.params
        self.params = self.form.dump()


    def attach(self):

        # Show the input form and save any changes that are made
        try:
            self.display_form()
        except Exception as e:
            st.exception(e)
            return

        # Catch any errors if incomplete inputs are provided
        if not self.form.complete:
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
    
    def save_changes_button(self):
        """
        When the "Save Changes" button is pressed,
        Save the values of the form to the MuData object,
        and also remove the 'edit-view' field from query_params.
        """

        # Save any changes which have been made to the form
        self.save_changes()

        # Clear the 'edit-view' query param
        del st.query_params['edit-view']

    def save_changes(self):
        """Save the values of the form to the MuData object."""

        # Get the views
        view = get_view(self.ix)

        # Set the params on the view
        view["params"] = self.form.dehydrate()

        set_view(self.ix, view)
