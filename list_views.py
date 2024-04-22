from mudata_explorer.helpers import all_views


for view in all_views:
    print("\n\t".join([
        view.name,
        view.type,
        view.desc,
        (
            "Categorie(s): " + ", ".join(view.categories)
            if view.categories
            else "No categories."
        )
    ]))
