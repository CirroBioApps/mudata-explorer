from cirro import DataPortalDataset
from cirro.sdk.file import DataPortalFile
from muon import MuData
import streamlit as st
from typing import Optional
from mudata_explorer.helpers.cirro_readers import util
from mudata_explorer.helpers.parsers import curatedMetagenomicData
from mudata_explorer.helpers.parsers import microbiome
from io import BytesIO

levels = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]


def read(
    dataset: DataPortalDataset
) -> Optional[MuData]:
    """Read datasets produced by curatedMetagenomicData."""

    # Select the file of interest
    file = select_file(dataset)

    # Parse the curatedMetagenomicData format as AnnData
    adata = curatedMetagenomicData.parse_df(file.read_csv(sep="\t"))

    # Run the microbiome analysis
    return microbiome.parse_adata(adata)


def select_file(dataset: DataPortalDataset) -> Optional[DataPortalFile]:
    # Get the list of files from the dataset
    files = util.list_files(dataset)

    # If there is more than one file, let the user select
    if len(files) > 1:
        selected_file = st.selectbox(
            "Select a file from the dataset",
            [f.name for f in files]
        )
        return next(
            file
            for file in files
            if file.name == selected_file
        )
    else:
        return files[0]
