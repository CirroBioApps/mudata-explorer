from mudata_explorer.base.view import View, ViewJSON
import streamlit as st


class Markdown(View):

    type = "markdown"
    name = "Markdown Text"
    help_text = "Write any text using markdown syntax."
    category = "Narrative"
    schema = {
        "text": {
            "type": "string",
            "label": "Text",
            "default": "",
            "help": "Use markdown syntax to write formatted text.",
            "multiline": True,
            "sidebar": True
        }
    }
    params = {}

    def display(self):

        text = self.params.get('text.value')
        if text is not None and len(text) > 0:
            st.write(text)

    def to_json(self) -> ViewJSON:
        return dict(
            type=self.type,
            figures=[self.params.get('text.value', '')]            
        )
