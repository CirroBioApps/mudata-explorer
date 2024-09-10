#!/usr/bin/env python3

from mudata_explorer.parsers import curatedMetagenomicData
from mudata_explorer.parsers import microbiome
from mudata_explorer.parsers.microbiome import MicrobiomeParams
from typing import Tuple
import json
from pathlib import Path
from pandas import read_csv, DataFrame
from anndata import AnnData
from mudata_explorer.sdk import io
import logging

logger = logging.getLogger("mudata-curated-metagenomic-data")


def run(adata: AnnData, params: MicrobiomeParams, basename: str):

    # Run the microbiome analysis
    mdata = microbiome._parse_adata(adata, groupby_var=True)
    microbiome._run_processes(mdata, params)
    microbiome._add_views(mdata, params)

    # Save the results
    io.write_h5mu(mdata, basename)


def pick_column(df) -> Tuple[str, bool]:
    """
    Select the column to compare by, and whether it is categorical
    """
    # Categorical data
    for cname in [
        'study_condition',
        'disease',
        'disease_subtype',
        'treatment',
        'non_westernized',
        'travel_destination',
        'body_subsite',
        'born_method',
        'anti_PD_1',
        'stec_count',
        'alcohol'
    ]:
        if has_multiple_groups(df, cname):
            return cname, True

    # Continuous data
    for cname in [
        "bmi",
        "age"
    ]:
        if has_multiple_groups(df, cname):
            return cname, False

    # Fallback categories
    for cname in [
        'visit_number',
        'age_category',
        'gender',
        'subject_id'
    ]:
        if has_multiple_groups(df, cname):
            return cname, True

    return None, None


def has_multiple_groups(df, cname):
    if cname in df.columns:
        vc = df[cname].value_counts()
        if (vc > 1).sum() > 1:
            return True
    return False


def setup_config(df: DataFrame, config: Path):
    """Set up a config file based on the contents of the dataset."""

    # Pick the metadata column to use to compare samples
    compare_by, is_categorical = pick_column(df)

    # Write the configuration file
    json.dump(
        [] if compare_by is None else [
            dict(
                dataset_name=(
                    config
                    .name
                    .replace(".config.json", "")
                    .replace("_", " ")
                ),
                compare_by=compare_by,
                label=compare_by.replace("_", " ").title(),
                n_top_features=20,
                is_categorical=is_categorical,
                leiden_res=1.0
            )
        ],
        config.open("w"),
        indent=4
    )


if __name__ == "__main__":

    # These studies do not have metadata annotations which facilitate
    # the comparison of the microbiome composition between groups

    for tsv in Path(".").rglob("*.relative_abundance.tsv"):

        basename = str(tsv).replace('.relative_abundance.tsv', '')

        done = Path(basename + '.done')
        if done.exists():
            continue

        # The config file contains information on how
        # the data should be processed
        config = Path(basename + '.config.json')

        # If there is no configuration file
        if not config.exists():

            # Read in the table
            df = read_csv(tsv, sep="\t")

            # If the dataset is at least 100 samples, then
            if df.shape[0] >= 100:

                # Set up a configuration file
                setup_config(df, config)

            # Otherwise, skip it
            else:
                continue

        # Read the configuration file
        config_list = [
            MicrobiomeParams(**cfg)
            for cfg in json.load(config.open())
        ]

        # Analyze each of the configured metadata categories
        for config_ix, config in enumerate(config_list):

            if config.compare_by is None:
                continue

            # Read the dataset as a DataFrame
            df = read_csv(tsv, sep="\t")

            # Drop any rows which are NaN for config.compare_by
            df = df.dropna(subset=[config.compare_by])

            # Convert to AnnData
            adata = curatedMetagenomicData.parse_df(df)

            # Name the file for the topic of analysis
            suffix = (
                config.label
                .replace(' ', '_')
                .replace('(', '')
                .replace(')', '')
                .lower()
            )

            # Run the analysis
            run(
                adata,
                config,
                # Name the output based on the column name
                f"{basename}-{config_ix}-{suffix}"
            )

        done.touch()
