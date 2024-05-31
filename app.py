#!/usr/bin/env streamlit run

from mudata_explorer import app
from mudata_explorer import pages

if __name__ == "__main__":

    # If no data has been uploaded
    if app.get_mdata() is None:

        # Let the user upload some data to get started
        pages.tables()

    # If there is data uploaded already
    else:

        # If there are figures to show
        if app.get_views() is not None:

            # Show the views
            pages.views()
