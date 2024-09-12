#!/usr/bin/env streamlit run

from mudata_explorer.app.mdata import list_modalities, get_mdata_exists
from mudata_explorer.helpers.views import get_views
from mudata_explorer import pages

if __name__ == "__main__":

    # If no data has been uploaded
    if (not get_mdata_exists()) or len(list_modalities()) == 0:

        # Let the user upload some data to get started
        pages.load()

    # If there is data uploaded already
    else:

        # If there are figures to show
        if len(get_views()) > 0:

            # Show the views
            pages.views()

        # If there are no figures to show
        else:

            # Show the tables
            pages.tables()
