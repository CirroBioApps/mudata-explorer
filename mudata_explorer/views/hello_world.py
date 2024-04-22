from mudata_explorer.base.view import View
from streamlit.delta_generator import DeltaGenerator


class HelloWorld(View):

    type = "hello-world"
    name = "Hello World"
    desc = "Basic hello world view."
    categories = ["Testing"]

    def display(self, container: DeltaGenerator):
        container.write("Hello, world!")
