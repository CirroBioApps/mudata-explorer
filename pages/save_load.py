from mudata_explorer import app
import streamlit as st
from streamlit.delta_generator import DeltaGenerator


def upload_button(container: DeltaGenerator):
    h5mu_file = container.file_uploader(
        "**Upload MuData (.h5mu)**"
    )
    if h5mu_file is None:
        return
    try:
        mdata = app.read_h5mu(h5mu_file)
    except ValueError as e:
        container.write(
            f"Could not parse file: {h5mu_file.name}\n\n{str(e)}"
        )
        return

    if app.hash_dat(app.mdata_to_binary(mdata)) in h5mu_file.name:
        container.write("**Data Validated**: Unique hash matches file name.")

    app.set_mdata(mdata)

    # Provide links to view the data, add more tables, or run processes
    app.show_shortcuts(
        [
            ("views", ":bar_chart: View Data"),
            ("processes", ":running: Run Processes"),
            ("add_data", ":page_facing_up: Add Tables")
        ],
        container=container
    )


def download_button(container: DeltaGenerator):

    # Get the current dataset along with its unique hash
    dat, hash = app.get_dat_hash()
    if dat is None:
        return

    # Compute the size of the file
    size = len(dat) / 1024
    # Format the size as a string
    if size < 1024:
        size = f"{size:.2f} KB"
    else:
        size = f"{size/1024:.2f} MB"

    container.write(f"**Unique Hash**: {hash}")
    container.write("---")

    # Name the downloaded file for the hash of the data
    container.write("**Download MuData**")
    container.download_button(
        f"Download {size}",
        dat,
        f"mudata-{hash}.h5mu",
        help="Click here to download the MuData object as a file."
    )


if __name__ == "__main__":

    app.setup_pages()

    container = st.container()
    container.write("#### Save / Load Dataset")
    container.write(
        "Download a snapshot of an existing dataset or upload a new one."
    )
    container.write("---")
    upload_button(container)
    container.write("---")
    download_button(container)
