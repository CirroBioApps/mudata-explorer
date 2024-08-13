from mudata_explorer.helpers.assets import all_views


for view in all_views:
    print("\n\t".join([
        view.name,
        view.type,
        (
            f"Category: {view.category}"
            if view.category is not None
            else "No categories."
        )
    ]))
