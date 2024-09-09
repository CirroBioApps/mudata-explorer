#!/usr/bin/env python3

from typing import Tuple
import json
from pathlib import Path
from pandas import read_csv, DataFrame
from muon import MuData
from mudata_explorer.sdk import io
from mudata_explorer.parsers.curatedMetagenomicData import read_tsv
from mudata_explorer.parsers.microbiome import MicrobiomeParams
from mudata_explorer.parsers.microbiome import _run_processes
from mudata_explorer.parsers.microbiome import _add_views
import logging
from scipy.spatial.distance import pdist
from scipy.stats import spearmanr
import umap

logger = logging.getLogger("mudata-curated-metagenomic-data")


def umap_params(df, metric="braycurtis", n_components=2):
    """
    Select the UMAP parameters which maximize the correlation
    of the embedded coordinates with the original data
    """

    # Get the pairwise distances for the input data
    true_dists = pdist(df, metric=metric)

    dat = []
    for min_dist in [0.1, 0.2, 0.3, 0.4]:
        for n_neighbors in [5, 10, 15, 20, 25]:
            # Run UMAP
            fit = umap.UMAP(
                metric=metric,
                n_neighbors=n_neighbors,
                min_dist=min_dist,
                n_components=n_components
            )
            umap_coords = fit.fit_transform(df)
            corr = spearmanr(
                true_dists,
                pdist(umap_coords, metric="euclidean")
            )
            dat.append(
                dict(
                    min_dist=min_dist,
                    n_neighbors=n_neighbors,
                    statistic=corr.statistic
                )
            )

    dat = DataFrame(dat).sort_values("statistic", ascending=False)

    return dict(
        metric=metric,
        n_components=n_components,
        min_dist=float(dat['min_dist'].values[0]),
        n_neighbors=int(dat['n_neighbors'].values[0])
    )


def leiden_params():
    return dict(
        metric="braycurtis",
        n_neighbors=15,
        resolution=1.,
    )


def run(mdata: MuData, config: MicrobiomeParams, basename: str):
    """
    1. Read in the TSV file as a DataFrame
    2. Separate out the metadata columns from the read counts
        (all of the read counts start with 'k__')
    3. Transform the read counts to proportions
    4. Build a MuData object with the metadata and transformed read counts
    5. Calculate summary metrics for the detected species
    6. Group the samples by microbiome composition
    7. Compare those groups to the groups defined in the disease column
    8. Find the organisms which are most differentially abundant between groups
    9. Group the samples by those organisms
    10. Compare those microbiome-based groups to the groups defined
    """

    # Run the analysis steps
    _run_processes(mdata, config)

    # Configure the visualizations
    _add_views(mdata, config)

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

        # Read the dataset
        mdata = read_tsv(tsv)

        # Read the configuration file
        config_list = [
            MicrobiomeParams(**cfg)
            for cfg in json.load(config.open())
        ]

        # Analyze each of the configured metadata categories
        for config_ix, config in enumerate(config_list):

            if config.compare_by is None:
                continue

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
                mdata,
                config,
                # Name the output based on the column name
                f"{basename}-{config_ix}-{suffix}"
            )

        done.touch()
