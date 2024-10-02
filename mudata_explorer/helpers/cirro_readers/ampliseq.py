from anndata import AnnData
from cirro import DataPortalDataset
from cirro.sdk.file import DataPortalFile
from mudata import MuData
import pandas as pd
import streamlit as st
from typing import Tuple, List, Optional
from mudata_explorer.helpers.cirro_readers import util
from mudata_explorer.parsers.microbiome import parse_adata

levels = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Species_exact"]


def read(
    dataset: DataPortalDataset
) -> Optional[MuData]:
    """Read datasets produced by nf-core/ampliseq."""

    adata = _read_ampliseq_as_anndata(dataset)
    if adata is None:
        return

    return parse_adata(adata, groupby_var=True)


def _read_ampliseq_as_anndata(dataset: DataPortalDataset) -> Optional[MuData]:
    # Get the list of files from the dataset
    files = util.list_files(dataset)

    # Read in the relative abundance table,
    # along with the taxonomic assignment for each feature
    with st.container(border=1):
        rel_abund = _read_rel_abund(files)
    if rel_abund is None:
        return
    tax, abund = rel_abund

    # Read the samplesheet saved with the dataset
    sample_meta = util.sample_metadata(dataset)

    # For this dataset, any dashes in sample names
    # are replaced with underscores
    sample_meta = sample_meta.rename(
        index=lambda i: str(i).replace("-", "_")
    )

    # Make an AnnData object
    adata = AnnData(
        X=abund.sort_index().T,
        var=tax.sort_index(),
        obs=sample_meta.reindex(abund.columns)
    )
    return adata


def _read_rel_abund(
    files: List[DataPortalFile]
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Read the relative abundance table from a dataset."""

    abund_file = util.find_file_by_extension(
        files,
        prefix="data/qiime2/rel_abundance_tables/rel-table-ASV_with-",
        suffix="-tax.tsv",
        selectbox_label="Use taxonomy assigned by:",
        none_found_msg="No relative abundance tables found."
    )
    if abund_file is None:
        return

    # Read in the table
    df = abund_file.read_csv(sep="\t", index_col=0)

    # Drop the unused columns
    df = (
        df
        .rename(columns=dict(confidence="Confidence"))
        .drop(columns=[
            col for col in ["Domain", "sequence", "Confidence"]
            if col in df.columns.values
        ])
    )

    # Parse the taxonomy
    tax = _parse_taxonomy(df)

    # Fix the species names
    tax = _fix_species(tax)

    # Remove any taxonomy inforamtion from the abundance table
    df = df.drop(columns=[l for l in levels + ["Taxon", "Confidence"] if l in df.columns.values])

    return tax, df


@st.cache_data
def _parse_taxonomy(df: pd.DataFrame):
    # If there is a "Taxon" column, parse it
    if "Taxon" in df.columns:

        df = pd.concat(
            [
                df.drop(columns=["Taxon"]),
                pd.DataFrame(
                    [
                        _parse_tax_string(taxon)
                        for taxon in df["Taxon"].values
                    ],
                    index=df.index
                )
            ],
            axis=1
        )

    # Only return the taxonomic ranks
    return df.reindex(columns=[l for l in levels if l in df.columns.values])


def _parse_tax_string(taxon: str):
    slc = dict(
        d="Kingdom",
        p="Phylum",
        c="Class",
        o="Order",
        f="Family",
        g="Genus",
        s="Species",
        t="Species_exact",
    )
    return {
        slc[part.strip()[0]]: part.strip()[3:]
        for part in taxon.split(";")
        if part != "Unassigned"
    }


def _fix_species(df: pd.DataFrame):
    for kw in ["Species", "Species_exact"]:
        if kw in df.columns:
            df = df.assign(**{
                kw: df.apply(
                    lambda r: (
                        f"{r['Genus']} {r[kw]}"
                        if r.get("Genus") and not pd.isnull(r[kw])
                        else r[kw]
                    ),
                    axis=1
                )
            })

    return df
