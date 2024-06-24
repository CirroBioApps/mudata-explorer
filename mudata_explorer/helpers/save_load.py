from mudata_explorer import app
from mudata_explorer.helpers import make_process
import streamlit as st
from streamlit.delta_generator import DeltaGenerator


@st.experimental_dialog("Upload MuData")
def upload_button():
    h5mu_file = st.file_uploader(
        "**Upload MuData (.h5mu)**"
    )
    if h5mu_file is None:
        return
    try:
        mdata = app.read_h5mu(h5mu_file)
    except ValueError as e:
        st.write(
            f"Could not parse file: {h5mu_file.name}\n\n{str(e)}"
        )
        return

    mdata_hash = app.hash_dat(h5mu_file.read())

    if mdata_hash in h5mu_file.name:
        st.write("**Data Validated**: Unique hash matches file name.")
    else:
        st.write("Unique file hash not found in file name.")

    if st.button("Load Dataset"):
        app.set_mdata(mdata)
        app.set_mdata_hash(mdata_hash)
        st.rerun()


def show_hash(container: DeltaGenerator):

    # Get the current dataset along with its unique hash
    dat, hash, _ = app.get_dat_hash()
    if dat is None:
        return

    container.write(f"**Unique Hash**: {hash}")


@st.experimental_dialog("Download MuData")
def download_button():

    # Get the current dataset along with its unique hash
    dat, hash, size = app.get_dat_hash()
    if dat is None:
        return

    st.write(f"File is {size}")
    st.write(f"Unique hash: {hash}")

    # Name the downloaded file for the hash of the data
    if st.download_button(
        "Download MuData (.h5mu)",
        dat,
        f"mudata-{hash}.h5mu",
        help="Click here to download the MuData object as a file."
    ):
        st.rerun()


@st.experimental_dialog("Select Previous Analysis", width="large")
def load_history():
    """Let the user select a previously run process to copy settings from."""

    # Get the list of previous analyses
    history = app.get_history(exclude=['add_data', 'add_view'])

    if len(history) == 0:
        st.write("No previous analyses found.")
        return

    # Instantiate a process for each one
    options = [
        make_process(h['process'], params=h['params'])
        for h in history
    ]

    # Make a display name for each
    display_names = [
        f"{process.name} ({h['timestamp']})"
        for process, h in zip(options, history)
    ]

    # Let the user select one of the previously-run processes
    selected_name = st.selectbox(
        "Select a previous analysis",
        display_names
    )
    if st.button("Load settings from selected analysis"):
        selected_ix = display_names.index(selected_name)
        proc = options[selected_ix]
        app.set_process(dict(
            type=proc.type,
            category=proc.category,
            params=proc.params
        ))
        st.rerun()
