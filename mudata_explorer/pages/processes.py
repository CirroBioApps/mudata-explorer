from typing import List
from mudata_explorer import app
from mudata_explorer.helpers import asset_categories, filter_by_category
from mudata_explorer.helpers import asset_type_desc_lists
from mudata_explorer.helpers import all_processes, make_process
import streamlit as st
from streamlit.delta_generator import DeltaGenerator


def update_process_kwargs(kw):
    return dict(
        on_change=app.update_process_on_change,
        key=f"process-{kw}",
        args=(kw,)
    )


def update_process_type(
    type_list: List[str],
    desc_list: List[str]
):
    selected_desc = st.session_state["process-type"]
    selected_process = type_list[desc_list.index(selected_desc)]
    app.update_process("type", selected_process)


def _get_selected_process(container: DeltaGenerator, process_def: dict):

    # Let the user select the type of process to add
    all_categories = asset_categories(all_processes)
    if process_def.get("category") not in all_categories:
        app.update_process("category", all_categories[0])

    selected_category = container.selectbox(
        "Select a category",
        all_categories,
        index=all_categories.index(process_def["category"]),
        **update_process_kwargs("category")
    )

    # Let the user select which process to add, filtering by category
    filtered_processes = filter_by_category(all_processes, selected_category)

    # Get the assets needed to select from the filtered views
    type_list, desc_list = asset_type_desc_lists(filtered_processes)

    if process_def.get("type") not in type_list:
        app.update_process("type", type_list[0])

    selected_desc = container.selectbox(
        "Select a process to run",
        desc_list,
        index=type_list.index(process_def["type"]),
        on_change=update_process_type,
        key="process-type",
        args=(type_list, desc_list,)
    )
    selected_process = type_list[desc_list.index(selected_desc)]
    return selected_process


def run():

    app.setup_sidebar(load_history=True)

    container = st.container()

    # Get the setup of the current process
    process_def = app.get_process()

    container.write("#### Analyze Data")

    if not app.has_mdata():
        container.page_link(
            "pages/tables.py",
            label="Upload data to get started"
        )
        return

    # Let the user select a process to run
    selected_process = _get_selected_process(container, process_def)

    # Instantiate the process using any parameters
    # which have been set
    process = make_process(
        selected_process,
        params=process_def.get("params", {})
    )

    # If there is any help text defined, show it
    if process.help_text:
        (
            container
            .expander("Description", expanded=False)
            .markdown(process.help_text)
        )

    # Get the parameters from the user
    process.get_data(container)

    # Get the output locations
    output_locs = process.get_output_locs()

    # If no data has been selected
    if len(output_locs) == 0:
        container.write("Please select input data")
        return

    # For each of the destination modalities
    for loc in output_locs:

        # Report to the user if data already exists in the destination key
        app.show_provenance(loc, container)

    # Let the user run the method, catching any errors
    if not process.params_complete:
        container.write("Please complete all input fields")
        return

    if container.button("Generate Results"):
        try:
            # Generate the results
            process.execute()
            st.rerun()

        except Exception as e:
            # Log the full traceback of the exception
            container.exception(e)
