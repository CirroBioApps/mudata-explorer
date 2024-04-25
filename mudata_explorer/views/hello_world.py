from mudata_explorer.base.view import View
from streamlit.delta_generator import DeltaGenerator


class HelloWorld(View):

    type = "hello-world"
    name = "Hello World"
    desc = "Basic hello world view."
    categories = ["Testing"]
    processed = True
    defaults = {"name": "World"}

    def display(self, container: DeltaGenerator):
        params = {
            kw: self.params.get(kw, val)
            for kw, val in self.defaults.items()
        }
        container.write(f"Hello, {params['name']}!")

    def inputs(self, form: DeltaGenerator):
        form.text_input(
            "Name",
            help="Enter your name here.",
            **self.param_kwargs("name")
        )
