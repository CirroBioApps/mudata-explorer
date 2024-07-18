import anndata as ad
from cirro import DataPortalDataset, DataPortalProject
from cirro.sdk.file import DataPortalFile
import streamlit as st
from streamlit.delta_generator import DeltaGenerator
from tempfile import TemporaryDirectory
from time import sleep
from typing import List, Optional
import re
import pandas as pd


def clear_cirro_client():
    # Wait for a second
    sleep(1)
    # Clear the client
    del st.session_state["Cirro"]
    # Let the user log in again
    st.rerun()


def list_files(
    dataset: DataPortalDataset,
    pattern: str = None
) -> List[DataPortalFile]:
    """List the files in a dataset."""

    # Try to get the list of files in the dataset
    try:
        files = dataset.list_files()
    # If there is an error
    except Exception as e:
        # Report it to the user
        st.error(f"Error: {e}")
        clear_cirro_client()

    # If a pattern was provided
    if pattern is not None:
        # Filter the files
        files = [f for f in files if re.match(pattern, f.name[5:])]

    return files


def select_file(
    container: DeltaGenerator,
    dataset: DataPortalDataset,
    pattern: str = None,
    files: Optional[List[DataPortalFile]] = None
) -> DataPortalFile:

    # If no files were provided
    if files is None:
        # Try to get the list of files in the dataset
        files = list_files(container, dataset)
    else:
        assert all(isinstance(f, DataPortalFile) for f in files)

    # Give the user the option to select a file
    # (strip off the 'data/' prefix)
    file = container.selectbox(
        "File",
        ["< select a file >"] + [
            f.name[5:] for f in files
            if pattern is None or re.match(pattern, f.name[5:])
        ],
        key="select-file"
    )

    if file == "< select a file >":
        return

    # Return the file object
    return next(f for f in files if f.name[5:] == file)


def sample_metadata(dataset: DataPortalDataset) -> List[dict]:
    # Get the project which the dataset is a part of
    project = DataPortalProject(
        dataset._client.projects.get(
            dataset.project_id
        ),
        dataset._client
    )
    # Return the sample metadata for that project
    return project.samples()


@st.cache_data(hash_funcs={DataPortalFile: lambda f: f.absolute_path})
def read_h5ad(h5ad_file: DataPortalFile) -> ad.AnnData:
    # Download the file to a temporary directory
    with TemporaryDirectory() as tmp:
        h5ad_file.download(tmp)
        filename = f"{tmp}/{h5ad_file.name}"
        return ad.read_h5ad(filename)


def guess_if_categorical(vals: pd.Series) -> bool:
    """
    Try to infer from a collection of values whether it is
    categorical or continuous in nature.
    True if categorical else False
    """
    vals = vals.dropna()
    if vals.apply(lambda v: isinstance(v, str)).any():
        return True

    vc = vals.value_counts()

    # If more than 25% of the values are singletons
    if (vc == 1).sum() > 0.25 * vals.shape[0]:
        return True

    return False
