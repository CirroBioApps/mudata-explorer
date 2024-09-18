from mudata_explorer.apps.microbiome import load
from mudata_explorer.apps.microbiome import menu

# Show the menu at the top
menu.show_menu(active_page="load_microbiome")

load.run()
