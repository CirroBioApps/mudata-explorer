import logging
from anndata import AnnData
import pandas as pd
from pathlib import Path
from mudata_explorer.sdk import io

logger = logging.getLogger(__name__)


def _parse_org_name(org_name: str) -> dict:
    """
    Parse the organism name to get the taxonomic ranks
    """

    ranks = {
        i[0]: i[3:]
        for i in org_name.split("|")
    }
    return {
        "kingdom": ranks['k'],
        "phylum": ranks['p'],
        "class": ranks['c'],
        "order": ranks['o'],
        "family": ranks['f'],
        "genus": ranks['g'],
        "species": ranks['s']
    }


def read_tsv(tsv: Path, sep="\t") -> AnnData:

    # Read in the DataFrame
    logger.info(f"Reading in {tsv}")
    df = pd.read_csv(tsv, sep=sep, index_col=0)

    return parse_df(df)


def parse_df(df: pd.DataFrame) -> AnnData:

    assert "sample_id" in df.columns, "No sample_id column found"
    df.set_index("sample_id", inplace=True)

    # Split up the metadata from the read counts
    obs = df.reindex(columns=[
        col for col in df.columns if not col.startswith("k__")
    ])
    counts = df.reindex(columns=[
        col for col in df.columns if col.startswith("k__")
    ])

    # Remove any columns from the metadata which are invariant
    obs = obs.reindex(
        columns=[
            cname for cname, cvals in obs.items()
            if cvals.nunique() > 1
        ]
    )

    # Get the taxonomic ranks from the species names
    tax_ranks = pd.DataFrame([
        _parse_org_name(org_name)
        for org_name in counts.columns
    ], index=counts.columns)

    # Strip down the species names to just the species
    counts = counts.rename(
        columns=lambda cname: cname.split("|")[-1][3:].replace("_", " ")
    )
    tax_ranks = tax_ranks.rename(
        index=lambda cname: cname.split("|")[-1][3:].replace("_", " ")
    )

    # Calculate the proportions from the counts
    abund = counts.apply(lambda r: r / r.sum(), axis=1)

    # Build the AnnData object
    adata = AnnData(abund, obs=obs, var=tax_ranks)

    adata.obs_names_make_unique()

    return adata
