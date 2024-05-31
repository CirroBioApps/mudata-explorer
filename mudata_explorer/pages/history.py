from typing import List
from mudata_explorer import app
import streamlit as st
from streamlit.delta_generator import DeltaGenerator


def print_history(history: List[dict], container: DeltaGenerator):
    for ix, event in enumerate(history):
        if "timestamp" in event:
            container.write(f"**{event['timestamp']}**")
        if "process" in event:
            container.write(f"- Ran: {event['process']}")
        if "updated_keys" in event:
            container.write(f"- Updated: {event['updated_keys']}")
        if "params" in event:
            for kw, val in event['params'].items():
                container.write(f" - {kw}: {val}")
        if ix < len(history) - 1:
            container.write("---")


def run():

    app.setup_sidebar()

    container = st.container()

    container.write("#### History")

    history = app.get_history()

    if len(history) == 0:
        container.write("No history to display.")

    else:
        print_history(history, container)
