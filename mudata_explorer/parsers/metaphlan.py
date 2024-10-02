from typing import Tuple
import pandas as pd


def parse_df(abund: pd.DataFrame, tax_level: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    # The samples are in the columns, so transpose the table
    abund = abund.T
    # Get the single letter code used to filter
    slc = tax_level[0].lower() if tax_level != "Strain" else "t"

    # Filter the table
    abund = abund[
        [
            col for col in abund.columns
            if col.split("|")[-1].startswith(f"{slc}__")
        ]
    ]

    # Make the full taxonomy table
    tax = pd.DataFrame(
        [
            {
                dict(
                    k="Kingdom",
                    p="Phylum",
                    c="Class",
                    o="Order",
                    f="Family",
                    g="Genus",
                    s="Species",
                    t="Strain"
                )[part[0]]: part.split("__", 1)[1].replace("_", " ")
                for part in col.split("|")
            }
            for col in abund.columns
        ],
        index=abund.columns
    )

    # Rename the taxa, both in the table and the taxonomy
    name_map = {
        col: [
            part.split("__")[1].replace("_", " ")
            for part in col.split("|")
            if part.startswith(f"{slc}__")
        ][0]
        for col in abund.columns
    }

    abund.rename(columns=name_map, inplace=True)
    tax.rename(index=name_map, inplace=True)
    
    # The values need to be converted to proportions
    abund = abund.apply(lambda x: x / x.sum(), axis=1)

    # Make sure that the number of dimensions adds up
    assert abund.shape[1] == tax.shape[0], "Number of taxa must match"

    return abund, tax
