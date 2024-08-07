from cirro import DataPortalDataset
from muon import MuData
from mudata_explorer import app
import streamlit as st
from tempfile import TemporaryDirectory
from typing import Optional
from mudata_explorer.helpers.cirro_readers import util
from muon import read_h5mu


def read(dataset: DataPortalDataset) -> Optional[MuData]:
    """Read from a dataset with 1 or more .h5mu files."""

    # Get the list of files from the dataset
    # filtering to only those with the .h5mu extension
    files = util.list_files(
        dataset,
        pattern=".*\.h5mu"
    )

    # If there is more than one file
    if len(files) > 1:

        # Pick the file from the dataset
        h5mu_file = util.select_file(
            dataset,
            files=files
        )
        # If no file is selected
        if h5mu_file is None:
            return

    # If there is only one file
    elif len(files) == 1:
        h5mu_file = files[0]
    # If there are no files
    else:
        st.error("No files found")
        return

    # Download to a temporary folder
    with TemporaryDirectory() as tmp:
        # Download the file
        try:
            h5mu_file.download(tmp)
        except Exception as e:
            st.error(f"Error downloading file: {e}")
            return
        # Path of the downloaded file
        filename = f"{tmp}/{h5mu_file.name}"
        # Read the MuData object
        mdata = read_h5mu(filename)
        # Calculate the hash of the object
        mdata_hash = app.hash_dat(open(filename, "rb").read())

    if mdata_hash in h5mu_file.name:
        st.write("**Data Validated**: Unique hash matches file name.")
    else:
        st.write("Unique file hash not found in file name.")

    return mdata
