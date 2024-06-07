from typing import List
from mudata_explorer import app
import streamlit as st
from streamlit.delta_generator import DeltaGenerator


def print_history(history: List[dict], container: DeltaGenerator):
    for event in history:
        label = f"**{event['timestamp']}**: Ran: {event['process']}"
        with container.expander(label, expanded=False):
            if "updated_keys" in event:
                st.write(f"- Updated: {event['updated_keys']}")
            if "params" in event:
                for kw, val in event['params'].items():
                    st.write(f" - {kw}: {val}")


def run():

    app.setup_sidebar()

    container = st.container()

    container.write("#### History")

    history = app.get_history()

    if len(history) == 0:
        container.write("No history to display.")

    else:
        print_history(history, container)
