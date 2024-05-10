from mudata_explorer import app
import streamlit as st

app.setup_pages()

settings = app.get_settings()

settings["editable"] = st.checkbox(
    "Edit Views",
    value=settings.get("editable", True),
    help="Allow the user to edit the views.",
)
app.show_shortcuts(
    [
        ("views", ":bar_chart: View Data")
    ]
)

app.set_settings(settings)
