from typing import List
from mudata_explorer import app
from mudata_explorer.helpers import plotting
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


def _get_selected_process(container: DeltaGenerator):

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


if __name__ == "__main__":

    app.setup_pages()

    container = st.container()

    settings = app.get_settings()

    # Get the setup of the current process
    process_def = app.get_process()

    container.write("#### Process Data")

    if app.has_mdata():

        # Let the user select a process to run
        selected_process = _get_selected_process(container)

        # Instantiate the process using any parameters
        # which have been set
        process = make_process(
            selected_process,
            params=process_def.get("params", {})
        )

        # Get the parameters from the user
        process.get_data(container)

        # Show any plots that the process generates
        process.display(container)

        container.write("##### Outputs")

        # Let the user select where to save the results

        # The user must select a modality
        all_modalities = app.list_modalities()
        dest_modality = container.selectbox(
            "Write results to modality:",
            all_modalities
        )

        # Get the key used for the resulting data
        dest_key = container.text_input(
            "Label to use for results:",
            process_def.get("dest_key", "results"),
            **update_process_kwargs("dest_key")
        )

        # Define the location where the output will be written
        loc = process.locate_results(
            dest_modality,
            dest_key
        )

        # Report to the user if data already exists in the destination key
        app.show_provenance(loc, container)

        # Let the user run the method, catching any errors
        if container.button("Generate Results"):
            try:

                # Generate the results
                res = process.execute()

                # Make sure that the results conform to the expected type
                assert isinstance(res, process.output_type), type(res)

                # Save the results
                process.save_results(loc, res)

            except Exception as e:
                # Log the full traceback of the exception
                container.exception(e)

            st.rerun()
