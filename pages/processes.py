from mudata_explorer import app
from mudata_explorer.helpers import asset_categories, filter_by_category
from mudata_explorer.helpers import asset_type_desc_lists
from mudata_explorer.helpers import all_processes, make_process
import streamlit as st


def update_process_kwargs(kw):
    return dict(
        on_change=app.update_process_on_change,
        key=f"process-{kw}",
        args=(kw,)
    )


if __name__ == "__main__":

    app.setup_pages()

    container = st.container()

    settings = app.get_settings()

    # Get the setup of the current process
    process = app.get_process()

    container.write("#### Process Data")

    # Let the user select the type of process to add
    all_categories = asset_categories(all_processes)
    if process.get("category") not in all_categories:
        app.update_process("category", all_categories[0])

    selected_category = container.selectbox(
        "Select a category",
        all_categories,
        index=all_categories.index(process["category"]),
        **update_process_kwargs("category")
    )

    # Let the user select which process to add, filtering by category
    filtered_processes = filter_by_category(all_processes, selected_category)

    # Get the assets needed to select from the filtered views
    type_list, desc_list = asset_type_desc_lists(filtered_processes)

    if process.get("type") not in type_list:
        app.update_process("type", type_list[0])

    selected_desc = container.selectbox(
        "Select a process to run",
        desc_list
    )
    selected_process = type_list[desc_list.index(selected_desc)]

    # Instantiate the process
    process = make_process(selected_process)

    # Run the process
    process.run(container)
