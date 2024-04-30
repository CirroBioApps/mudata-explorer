from mudata_explorer import app
import streamlit as st

app.setup_pages()

st.write("#### Settings")

settings = app.get_settings()

st.write("---")
st.write("#### Editable Views")
settings["editable"] = st.checkbox(
    "Enable",
    value=settings.get("editable", True),
    help="Allow the user to edit the views.",
)

app.set_settings(settings)
