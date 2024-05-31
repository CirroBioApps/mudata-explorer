from mudata_explorer import app
import streamlit as st


def run():
    app.setup_sidebar()

    if app.get_mdata() is None:
        st.page_link(
            "pages/tables.py",
            label="Upload data to get started"
        )
        return

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
