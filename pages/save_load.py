import json
import muon as mu
from mudata_explorer import app
import streamlit as st
from streamlit.delta_generator import DeltaGenerator
from tempfile import NamedTemporaryFile


def upload_button(container: DeltaGenerator):
    h5mu_file = container.file_uploader(
        "**Upload MuData (.h5mu)**"
    )
    if h5mu_file is None:
        return
    try:
        mdata = read_h5mu(h5mu_file)
    except ValueError as e:
        container.write(
            f"Could not parse file: {h5mu_file.name}\n\n{str(e)}"
        )
        return

    app.set_mdata(mdata)


def read_h5mu(h5mu_file):
    with open(h5mu_file.name, "wb") as f:
        f.write(h5mu_file.getvalue())
    return mu.read(h5mu_file.name)


def download_button(container: DeltaGenerator):

    app.summarize_mdata(container)

    mdata = app.get_mdata()
    if mdata is None:
        return

    # Convert the view data to a JSON string
    mdata.uns["mudata-explorer-views"] = json.dumps(
        mdata.uns["mudata-explorer-views"]
    )

    # Write out the MuData object to a temporary file
    with NamedTemporaryFile(suffix=".h5mu", delete=True) as tmp:
        # Format as h5mu
        mdata.write(tmp.file.name)

        # Get the file object in bytes
        dat = tmp.read()

        # Compute the hash of the file
        hash = app.hash_dat(dat)

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
    container.write("### Save / Load Dataset")
    container.write("Download a snapshot of an existing dataset or upload a new one.")
    container.write("---")
    upload_button(container)
    container.write("---")
    download_button(container)
