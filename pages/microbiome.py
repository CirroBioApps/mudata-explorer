from mudata_explorer.apps.microbiome import menu
from mudata_explorer.apps.microbiome import load_explanation
from mudata_explorer.app.view import view_non_editable

# Load the explanatory dataset
load_explanation()

# Show the menu at the top
menu.show_menu(active_page="microbiome")

# Display the explanation
view_non_editable()
