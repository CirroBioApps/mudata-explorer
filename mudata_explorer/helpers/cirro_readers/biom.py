from tempfile import TemporaryDirectory
from biom import load_table
from anndata import AnnData
from cirro import DataPortalDataset
from mudata import MuData
import streamlit as st
from typing import Optional
from mudata_explorer.helpers.cirro_readers import util
from mudata_explorer.parsers.microbiome import parse_adata


def read(
    dataset: DataPortalDataset
) -> Optional[MuData]:
    """Read BIOM datasets."""

    adata = _read_biom_as_anndata(dataset)
    if adata is None:
        return
    
    return parse_adata(adata)


def _read_biom_as_anndata(
    dataset: DataPortalDataset
) -> Optional[AnnData]:
    """Read the BIOM file from a dataset."""

    # Get the list of files from the dataset
    files = util.list_files(dataset)

    biom_file = util.find_file_by_extension(
        files,
        prefix="data/",
        suffix=".biom",
        selectbox_label="Use BIOM file:",
        none_found_msg="No BIOM files found."
    )
    if biom_file is None:
        return

    # Download the file and read it
    with st.spinner("Reading BIOM file"):
        with TemporaryDirectory() as tmpdir:
            biom_file.download(tmpdir)
            table = load_table(f"{tmpdir}/{biom_file.name}")
    return table.to_anndata()
