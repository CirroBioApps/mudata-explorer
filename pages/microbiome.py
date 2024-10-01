from mudata_explorer.apps.microbiome import menu, explanation
from mudata_explorer.app.view import view_non_editable
from mudata_explorer.app.mdata import set_mdata
from copy import copy

# Load the explanatory dataset
set_mdata(copy(explanation), full=True)

# Show the menu at the top
menu.show_menu(active_page="microbiome")

# Display the explanation
view_non_editable()
