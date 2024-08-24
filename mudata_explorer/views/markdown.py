from mudata_explorer.base.view import View
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

    def display(self):

        text = self.params.get('text')
        if text is not None and len(text) > 0:
            st.write(text)
