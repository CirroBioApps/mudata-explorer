import json
from mudata_explorer import app
import streamlit as st


def print_history(history, container):
    container.write(json.dumps(history, indent=4))


if __name__ == "__main__":

    app.setup_pages()

    container = st.container()

    container.write("#### History")

    history = app.get_history()

    if len(history) == 0:
        container.write("No history to display.")

    else:
        print_history(history, container)
