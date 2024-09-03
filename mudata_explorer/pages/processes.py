import pandas as pd
from mudata_explorer.helpers.assets import asset_categories, filter_by_category
from mudata_explorer.helpers.assets import asset_dataframe
from mudata_explorer.helpers.assets import all_processes, make_process
from mudata_explorer.app.process import get_process, update_process, update_process_on_change
from mudata_explorer.app.sidebar import setup_sidebar
from mudata_explorer.app.provenance import show_provenance
from mudata_explorer.app.mdata import has_mdata
import streamlit as st
from streamlit.delta_generator import DeltaGenerator


def update_process_kwargs(kw):
    return dict(
        on_change=update_process_on_change,
        key=f"process-{kw}",
        args=(kw,)
    )


def update_process_type(df: pd.DataFrame):
    selected_name = st.session_state["process-name"]
    selected_process = df["type"].values[
        df["name"].tolist().index(selected_name)
    ]
    update_process("type", selected_process)


def _get_selected_process(container: DeltaGenerator, process_def: dict):

    # Let the user select the type of process to add
    all_categories = asset_categories(all_processes)
    if process_def.get("category") not in all_categories:
        update_process("category", all_categories[0])

    selected_category = container.selectbox(
        "Select a category",
        all_categories,
        index=all_categories.index(process_def["category"]),
        **update_process_kwargs("category")
    )

    # Let the user select which process to add, filtering by category
    filtered_processes = filter_by_category(all_processes, selected_category)

    # Get the assets needed to select from the filtered views
    df = asset_dataframe(filtered_processes)

    if process_def.get("type") not in df["type"].values:
        update_process("type", df["type"].values[0])

    selected_name = container.selectbox(
        "Select a process to run",
        df["name"].values,
        index=df["type"].tolist().index(process_def["type"]),
        on_change=update_process_type,
        key="process-name",
        args=(df,)
    )
    selected_ix = df["name"].tolist().index(selected_name)
    selected_process = df["type"].values[selected_ix]
    return selected_process


def run():

    setup_sidebar(load_history=True, page_layout="wide")

    container = st.container()

    # Get the setup of the current process
    process_def = get_process()

    container.write("#### Analyze Data")

    if not has_mdata():
        container.page_link(
            "pages/tables.py",
            label="Upload data to get started",
            icon=":material/table:"
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
    with container:
        process.form.render()
    process.params = process.form.dump()

    # Get the output locations
    output_locs = process.get_output_locs()

    # If no data has been selected
    if len(output_locs) == 0:
        container.write("Please select input data")
        return

    # Check if the data already exists in the destination key
    if any([
        # Report to the user if data already exists in the destination key
        show_provenance(loc, container)
        # For each of the destination modalities
        for loc in output_locs
    ]):
        container.write(
            "> Existing data will be replaced"
        )

    # Let the user run the method, catching any errors
    if not process.form.complete:
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
