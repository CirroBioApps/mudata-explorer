from anndata import AnnData
from cirro import DataPortalDataset
from cirro.sdk.file import DataPortalFile
from mudata import MuData
import pandas as pd
import streamlit as st
from typing import Tuple, List, Optional
from mudata_explorer.helpers.cirro_readers import util
from mudata_explorer.parsers import metaphlan
from mudata_explorer.parsers.microbiome import parse_adata


def read(
    dataset: DataPortalDataset
) -> Optional[MuData]:
    """Read datasets produced by nf-core/ampliseq."""

    adata = _read_metaphlan_as_anndata(dataset)
    if adata is None:
        return

    return parse_adata(adata, groupby_var=False)


def _read_metaphlan_as_anndata(dataset: DataPortalDataset) -> Optional[AnnData]:
    # Get the list of files from the dataset
    files = util.list_files(dataset)

    # Read in the relative abundance table,
    # which also contains the taxonomic assignment for each feature
    with st.container(border=1):
        rel_abund = _read_rel_abund(files)
    if rel_abund is None:
        return
    abund, tax = rel_abund

    # Read the samplesheet saved with the dataset
    sample_meta = util.sample_metadata(dataset)

    # Make an AnnData object
    adata = AnnData(
        X=abund,
        var=tax,
        obs=sample_meta
    )
    return adata


def _read_rel_abund(
    files: List[DataPortalFile]
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Read the relative abundance table from a dataset."""

    abund_file = util.find_file_by_extension(
        files,
        prefix="data/",
        suffix="merged_abundance_table.tsv",
        selectbox_label="Use taxonomy assigned by:",
        none_found_msg="No file found: merged_abundance_table.tsv"
    )
    if abund_file is None:
        return

    # Read in the table
    abund = abund_file.read_csv(sep="\t", index_col=0, comment="#")

    # The user must select which taxonomic level to use
    tax_level = st.selectbox(
        "Taxonomic Level",
        ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"],
        index=6
    )

    # Parse the table using the syntax of MetaPhlAn
    return metaphlan.parse_df(abund, tax_level)
