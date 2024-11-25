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
    """Read datasets produced by gig-map/align_pangenome.nf."""

    adata = _read_gig_map_as_anndata(dataset)
    if adata is None:
        return

    return parse_adata(adata, groupby_var=False, sum_to_one=False)


def _read_gig_map_as_anndata(dataset: DataPortalDataset) -> Optional[AnnData]:
    # Get the list of files from the dataset
    files = util.list_files(dataset)

    # Read in the relative abundance table,
    # which also contains the taxonomic assignment for each feature
    with st.container(border=1):
        feature_type = {
            "Gene Bins": "bins",
            "Genes": "genes",
        }[
            st.selectbox(
                "Feature Type",
                ["Gene Bins", "Genes"],
            )
        ]
        filename = f"metagenome.{feature_type}.h5ad"

        abund_file = util.find_file_by_extension(
            files,
            prefix="data/",
            suffix=filename,
            selectbox_label="Use data from:",
            none_found_msg=f"No file found: {filename}"
        )
        if abund_file is None:
            return

        # Read in the table
        abund = abund_file.read_h5ad()

        # Remove the NaN feature (genes which are not in any bin)
        abund = abund[:, [x for x in abund.var.index if not x == "nan"]]

        # If the gene bins were selected
        if feature_type == "bins":
            # Optionally weight all abundance values by the number of genes in each bin
            if st.checkbox("Weight abundances by the number of genes in each bin", value=True):
                abund.X = (abund.X / abund.var["n_genes"].values) * abund.var["n_genes"].sum()

        # Ask the user how they want the sequencing depth to be normalized
        # Options:
        # - Unadjusted (raw counts)
        # - Proportion of Total Reads
        # - Proportion of Aligned Reads
        # - Relative to Selected Bin (or Gene)

        norm_method = st.selectbox(
            "Sequencing Depth Normalization Method:",
            [
                "Proportion of Total Reads",
                "Proportion of Aligned Reads",
                "Relative to Selected Bin (or Gene)",
                "None (raw counts)",
            ]
        )
        if norm_method == "Proportion of Total Reads":
            abund.X = (abund.X.T / abund.obs["genes:tot_reads"].values).T
            st.write("Sum of values for each sample is <= 1")
        elif norm_method == "Proportion of Aligned Reads":
            abund.X = (abund.X.T / abund.obs["genes:aligned_reads"].values).T
            st.write("Sum of values for each sample is ~= 1")
        elif norm_method == "Relative to Selected Bin (or Gene)":
            # Ask the user to select a bin or gene
            selected_feature = st.selectbox(
                "Select a Bin or Gene:",
                abund.var.index
            )
            abund.X = (abund.X.T / abund.to_df()[selected_feature].values).T
            # Remove the selected feature from the table
            abund = abund[:, [x for x in abund.var.index if not x == selected_feature]]
            st.write("Sum of values for each sample is > 1 or < 1")
        else:
            st.write("Sum of values for each sample is the number of aligned reads")

        return abund
