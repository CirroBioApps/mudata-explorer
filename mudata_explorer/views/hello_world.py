from mudata_explorer.base.view import View
import streamlit as st


class HelloWorld(View):

    type = "hello-world"
    name = "Hello World"
    help_text = "Basic hello world view."
    category = "Testing"
    defaults = {"name": "World"}

    def display(self):
        if self.params_editable:
            self.params["name"] = st.text_input(
                "Name",
                help="Enter your name here.",
                **self.input_value_kwargs("name")
            )

        st.write(f"Hello, {self.params['name']}!")
