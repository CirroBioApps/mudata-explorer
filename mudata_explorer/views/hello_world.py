from mudata_explorer.base.view import View
from streamlit.delta_generator import DeltaGenerator


class HelloWorld(View):

    type = "hello-world"
    name = "Hello World"
    help_text = "Basic hello world view."
    category = "Testing"
    defaults = {"name": "World"}

    def display(self, container: DeltaGenerator):
        if self.params_editable:
            self.params["name"] = container.text_input(
                "Name",
                help="Enter your name here.",
                **self.input_value_kwargs("name")
            )

        container.write(f"Hello, {self.params['name']}!")
