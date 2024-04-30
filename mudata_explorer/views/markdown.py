from mudata_explorer.base.view import View
from streamlit.delta_generator import DeltaGenerator


class Markdown(View):

    type = "markdown"
    name = "Markdown Text"
    desc = "Write any text using markdown syntax."
    categories = ["Narrative"]
    defaults = {"text": ""}

    def display(self, container: DeltaGenerator):
        if self.editable:
            self.params["text"] = self.inputs_container.text_area(
                "Markdown Text",
                help="Use markdown syntax to write formatted text.",
                **self.input_value_kwargs("text")
            )

        container.write(self.params['text'])
