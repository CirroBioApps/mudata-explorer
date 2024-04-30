from mudata_explorer import app
import streamlit as st

if __name__ == "__main__":

    app.setup_pages()

    container = st.container()

    container.write("#### History")

    history = app.get_history()

    if len(history) == 0:
        container.write("No history to display.")

    else:
        container.write("\n\n".join(history))
