#!/usr/bin/env streamlit run

from mudata_explorer import app
from mudata_explorer import pages

if __name__ == "__main__":

    # If no data has been uploaded
    if app.get_mdata() is None:

        # Let the user upload some data to get started
        pages.tables()

    # If no observations have been uploaded
    elif len(app.list_modalities()) == 0:

        # Let the user upload some data to get started
        pages.tables()

    # If there is data uploaded already
    else:

        # If there are figures to show
        if isinstance(app.get_views(), list) and len(app.get_views()) > 0:

            # Show the views
            pages.views()

        # If there are no figures to show
        else:

            # Show the processes
            pages.processes()
