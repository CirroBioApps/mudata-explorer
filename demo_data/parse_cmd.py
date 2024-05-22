#!/usr/bin/env python3

import click
import pandas as pd


@click.command()
@click.argument('input_file', type=click.Path(exists=True))
def read_tsv(input_file):

    assert input_file.endswith(".tsv")
    # Read the TSV file
    df = pd.read_csv(input_file, sep='\t', index_col=0)

    output_prefix = input_file.split("/")[-1].replace(".tsv", "")

    # Split up the metadata and the data
    metadata_columns = [
        cname for cname in df.columns.values
        if not cname[1:3] == '__'
    ]
    (
        df[metadata_columns]
        .to_csv(
            output_prefix + "_metadata.csv"
        )
    )
    (
        df
        .drop(columns=metadata_columns)
        .to_csv(
            output_prefix + "_data.csv"
        )
    )


if __name__ == "__main__":
    read_tsv()
