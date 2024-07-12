#!/usr/bin/env python3

from typing import Optional
from dataclasses import dataclass
import json
from pathlib import Path
from pandas import read_csv, DataFrame, qcut
from muon import MuData
from mudata_explorer.sdk import view, io, process
import logging
from scipy.spatial.distance import pdist
from scipy.stats import spearmanr
import umap

logger = logging.getLogger("mudata-curated-metagenomic-data")


@dataclass
class Config:
    """Configuration for the analysis."""
    # Name of the dataset, to be used in the title
    dataset_name: str
    # Column name to use for comparison
    cname: str
    # Label to use when describing those comparison groups
    label: str
    # Number of features to use in the analysis and plots
    n_features: int
    # Optional transformation to apply to the column
    transform: Optional[str] = None


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


def transform(mdata: MuData, config: Config):
    if config.transform == "quartiles":
        kw = f"{config.cname}_quartiles"
        mdata.obs[kw] = qcut(
            mdata.obs[config.cname],
            q=4,
            labels=[
                f"Q{i+1}"
                for i in range(4)
            ]
        )
        config.cname = kw


def run(mdata: MuData, config: Config, basename: str):
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

    # Transform the metadata column if necessary
    transform(mdata, config)

    # Run the analysis steps
    analysis(mdata, config)

    # Configure the visualizations
    views(mdata, config)

    # Save the results
    io.write_h5mu(mdata, basename)


def read_mdata(tsv: Path) -> MuData:

    # Read in the DataFrame
    logger.info(f"Reading in {tsv}")
    df = read_csv(tsv, sep="\t", index_col=0)

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
    tax_ranks = DataFrame([
        parse_org_name(org_name)
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
    prop = counts / counts.sum()

    # Build the MuData object
    mdata = io.build_mdata(
        dict(
            # counts=counts,
            prop=prop
        ),
        obs=obs
    )

    mdata.mod["prop"].var = tax_ranks
    return mdata


def parse_org_name(org_name: str) -> dict:
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


def analysis(mdata: MuData, config: Config):
    """
    Run the analysis steps on the MuData object
    """

    # Steps 5-10
    summary_metrics(mdata)
    group_samples_unweighted(mdata)
    diff_abund(mdata, "leiden")
    ordinate_unweighted(mdata)

    diff_abund(mdata, config.cname)
    group_samples_weighted(mdata, config.cname, n_features=config.n_features)
    diff_abund(mdata, f"leiden_{config.cname}")
    ordinate_weighted(mdata, config.cname, n_features=config.n_features)

    # Run UMAP on the features (species)
    ordinate_features(mdata)


def views(mdata: MuData, config: Config):
    """
    1. Add a header markdown element
    2. Show the unweighted UMAP, colored by the Leiden groups
        and then also by the metadata column of interest
    """

    header(mdata, config.dataset_name)
    text(mdata, "### Differentially abundant species across community types")

    scatter_umap_unweighted(mdata, "leiden", label="Cluster (Unweighted)")
    text(mdata, """**Figure 1. UMAP ordination of samples colored by community type.** 
The overall microbiome composition of the samples is represented by UMAP ordination,
where the samples are colored by the groups (i.e. 'community types') which were
determined by unsupervised clustering (Leiden).""") # noqa

    dotplot_assoc_features(mdata, "leiden", config.n_features, label="Cluster")
    text(mdata, f"""**Figure 2. Relative abundance of organisms which vary
across community types.**
The relative abundance of the top {config.n_features:,} species which
are most differentially abundant between those community types
is shown using a dotplot display.
The size of each point represents the average abundance of
that species across all samples from each group of samples
indicated on the y-axis.
To account for different scales of abundance, the Z-score normalization
is used in which the mean is subtracted from each value and then divided
by the standard deviation.""")

    text(mdata, f"### Microbiome Composition ~ {config.label}")

    scatter_umap_unweighted(mdata, config.cname, config.label)
    text(mdata, f"""**Figure 3. UMAP ordination of samples colored by {config.label}.**
The overall microbiome composition of the samples is represented by UMAP ordination,
where the samples are colored by the '{config.label}' metadata category.""") # noqa

    compare_groups(mdata, "leiden", "Cluster (Unweighted)", config.cname, config.label)
    text(mdata, f"""**Figure 4. Comparison of community types to {config.label} groups.**
The barplot shows the number of samples in each Leiden group, broken up by
which are also in each '{config.label}' group.""") # noqa

    text(mdata, f"### Microbes Associated with {config.label} Groups")

    scatter_umap_weighted(mdata, config.cname, config.cname, config.label)
    text(mdata, f"""**Figure 5. UMAP ordination of samples weighted by {config.label}-associated microbes.**
Only the {config.n_features:,} species which are most differentially
abundant between the {config.label} groups are used to ordinate the samples.
""") # noqa

    dotplot_assoc_features(
        mdata,
        config.cname,
        config.n_features,
        label=config.label
    )
    text(mdata, f"""**Figure 6. Relative abundance of organisms which vary across {config.label} groups.**
The dotplot shows the normalized abundance of the top {config.n_features} species which
are most differentially abundant between the {config.label} groups.""") # noqa

    scatter_assoc_features(mdata, config)
    text(mdata, f"""**Figure 7. UMAP projection of organisms which vary across {config.label} groups.**
The scatterplot shows the UMAP projection of the samples colored by the taxonomic class.
The distance between each point reflects the similarity of the microbial composition
across samples from the complete dataset.
Only the top {config.n_features} species which are most differentially abundant between the
{config.label} groups are shown.
The size of each point represents the F-Statistic, indicating the relative degree of
difference in abundance between groups of samples.""") # noqa

    text(mdata, f"### Community Types Associated with {config.label} Groups")
    scatter_umap_weighted(
        mdata,
        config.cname,
        f"leiden_{config.cname}",
        f"Cluster ({config.label}-Weighted)"
    )
    text(mdata, f"""**Figure 8. Community types weighted by {config.label}-associated microbes.**
Only the {config.n_features:,} species which are most differentially
abundant between the {config.label} groups are used to group the samples into community types
using Leiden clustering.""") # noqa

    dotplot_assoc_features(
        mdata,
        f"leiden_{config.cname}",
        config.n_features,
        label=f"Cluster ({config.label}-weighted)"
    )
    text(mdata, f"""**Figure 9. Relative abundance of organisms which vary
across {config.label} groups.** The dotplot shows the normalized abundance of
the top {config.n_features} species which are most differentially abundant
between the {config.label} groups.""")

    compare_groups(
        mdata,
        f"leiden_{config.cname}",
        f"Cluster ({config.label}-weighted)",
        config.cname,
        config.label
    )
    text(mdata, f"""**Figure 10. Comparison of {config.label}-weighted
community types to {config.label} groups.** The barplot shows the
number of samples in each Leiden group, broken up by which are also
in each '{config.label}' group.""")


def header(mdata, name):

    text(
        mdata,
        f"""# {name}
### Curated Metagenomic Data"""
    )


def text(mdata, text):
    view.markdown(mdata, text=text)


def scatter_umap_unweighted(mdata: MuData, color_by: str, label: str):
    view.plotly_scatter(
        mdata,
        data_axis=0,
        data_color=dict(
            table=["Observation Metadata"],
            cname=color_by,
            label=label,
            is_categorical=True,
            scale="D3",
            enabled=True
        ),
        data_size_enabled=False,
        **{
            f"data_{ax}": dict(
                table=["prop.obsm.umap"],
                label=f"UMAP {n}",
                cname=f"UMAP {n}"
            )
            for ax, n in zip(["x", "y"], [1, 2])
        }
    )


def scatter_umap_weighted(
    mdata: MuData,
    cname: str,
    color_by: str,
    color_label: str
):
    view.plotly_scatter(
        mdata,
        data_axis=0,
        data_color=dict(
            table=["Observation Metadata"],
            cname=color_by,
            label=color_label,
            is_categorical=True,
            scale="D3",
            enabled=True
        ),
        data_size_enabled=False,
        **{
            f"data_{ax}": dict(
                table=[f"prop.obsm.umap_{cname}"],
                label=f"UMAP {n}",
                cname=f"UMAP {n}"
            )
            for ax, n in zip(["x", "y"], [1, 2])
        }
    )


def scatter_assoc_features(
    mdata: MuData,
    config: Config,
):
    table = f"prop.varm.kruskal_{config.cname}"
    view.plotly_scatter(
        mdata,
        data=dict(
            axis=1,
            color=dict(
                table=["prop.metadata"],
                cname="class",
                label="Taxonomic Class",
                is_categorical=True,
                scale="D3",
                enabled=True
            ),
            size=dict(
                table=[table],
                cname="f_statistic",
                label=f"{config.label} Association (F-Statistic)",
                enabled=True
            ),
            rows_query_query=dict(
                table=[table],
                cname="rank",
                type="value",
                expr="<=",
                value=str(config.n_features)
            ),
            **{
                ax: dict(
                    table=["prop.varm.umap"],
                    label=f"UMAP {n}",
                    cname=f"UMAP {n}"
                )
                for ax, n in zip(["x", "y"], [1, 2])
            }
        )
    )


def dotplot_assoc_features(mdata: MuData, cname: str, n_features=20, label=None):
    view.plotly_category_summarize_values(
        mdata,
        **{
            "table": {
                "category": {
                    "category": {
                        "table": [
                            "Observation Metadata"
                        ],
                        "cname": cname,
                        "label": label if label else cname
                    }
                },
                "data": {
                    "tables": [
                        "prop.data"
                    ],
                    "cols_query": {
                        "query": {
                            "cname": "rank",
                            "type": "value",
                            "table": [
                                f"prop.varm.kruskal_{cname}"
                            ],
                            "expr": "<=",
                            "value": str(n_features)
                        }
                    },
                    "transforms": [
                        "zscores_cols"
                    ]
                }
            },
            "formatting": {
                "color": "None",
                "sort_by": "Values"
            }
        }
    )


def boxplot_assoc_features(mdata: MuData, cname: str, n_features=20):
    view.plotly_box_multiple(
        mdata,
        **{
            "table": {
                "category": {
                    "category": {
                        "table": [
                            "Observation Metadata"
                        ],
                        "cname": cname,
                        "label": cname
                    }
                },
                "data": {
                    "tables": [
                        "prop.data"
                    ],
                    "cols_query": {
                        "query": {
                            "cname": "rank",
                            "type": "value",
                            "table": [
                                f"prop.varm.kruskal_{cname}"
                            ],
                            "expr": "<=",
                            "value": str(n_features)
                        }
                    }
                }
            },
            "variable_options": {
                "sort_by": "Name",
                "axis": "Y-Axis"
            },
            "category_options": {
                "axis": "Color"
            },
            "display_options": {
                "var_label": "",
                "val_label": "Relative Abundance"
            }
        }
    )


def compare_groups(
    mdata: MuData,
    group1: str,
    label1: str,
    group2: str,
    label2: str
):
    """
    Make a barplot comparing the groups in group1 to the groups in group2
    """
    view.plotly_category_count(
        mdata,
        **{
            "barmode": "group",
            "data": {
                "color": {
                    "cname": group2,
                    "label": label2,
                    "table": [
                        "Observation Metadata"
                    ]
                },
                "x": {
                    "cname": group1,
                    "label": label1,
                    "table": [
                        "Observation Metadata"
                    ]
                }
            },
            "annotation_options.chisquare": True
        }
    )


def pick_column(df) -> str:
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
            return cname, None

    # Continuous data
    for cname in [
        "bmi",
        "age"
    ]:
        if has_multiple_groups(df, cname):
            return cname, "quartiles"

    # Fallback categories
    for cname in [
        'visit_number',
        'age_category',
        'gender',
        'subject_id'
    ]:
        if has_multiple_groups(df, cname):
            return cname, None

    return None, None


def has_multiple_groups(df, cname):
    if cname in df.columns:
        vc = df[cname].value_counts()
        if (vc > 1).sum() > 1:
            return True
    return False


def summary_metrics(mdata: MuData):
    """
    Calculate summary metrics for the detected species
    """

    process.summary_stats(
        mdata,
        outputs_dest_key="summary_stats",
        table_data_axis=1,
        table_data_tables=["prop.data"]
    )


def group_samples_unweighted(
    mdata: MuData,
):
    group_samples(
        mdata,
        outputs_dest_key="leiden",
        **filter_prevalence("prop")
    )


def group_samples_weighted(
    mdata: MuData,
    metadata_category: str,
    n_features=20
):
    group_samples(
        mdata,
        outputs_dest_key=f"leiden_{metadata_category}",
        **filter_diff_abund("prop", metadata_category, n_features)
    )


def filter_prevalence(mod: str, min_prop=0.1):
    return dict(
        table_data_axis=0,
        table_data_tables=[f"{mod}.data"],
        table_data_cols_query_query_table=[f"{mod}.varm.summary_stats"],
        table_data_cols_query_query_type="value",
        table_data_cols_query_query_cname="prop_positive",
        table_data_cols_query_query_expr=">=",
        table_data_cols_query_query_value=str(min_prop),
    )


def filter_prevalence_df(mdata: MuData, mod: str, min_prop=0.1):
    adata = mdata.mod[mod]

    return adata[
        :,
        adata.varm['summary_stats']['prop_positive'] >= min_prop
    ].to_df()


def filter_diff_abund(
    mod: str,
    metadata_category: str,
    n_features: int
):
    return dict(
        table_data_axis=0,
        table_data_tables=[f"{mod}.data"],
        table_data_cols_query_query_table=[
            f"{mod}.varm.kruskal_{metadata_category}"
        ],
        table_data_cols_query_query_type="value",
        table_data_cols_query_query_cname="rank",
        table_data_cols_query_query_expr="<",
        table_data_cols_query_query_value=str(n_features)
    )


def filter_diff_abund_df(
    mdata: MuData,
    mod: str,
    metadata_category: str,
    n_features: int
):
    adata = mdata.mod[mod]

    return adata[
        :,
        adata.varm[f'kruskal_{metadata_category}']['rank'] <= n_features
    ].to_df()


def group_samples(
    mdata: MuData,
    outputs_dest_key="leiden",
    **kwargs
):
    """
    Group the samples by microbiome composition
    """

    # Run Leiden clustering
    process.leiden(
        mdata,
        clustering=leiden_params(),
        outputs_dest_key=outputs_dest_key,
        **kwargs
    )


def ordinate_unweighted(mdata: MuData):
    ordinate(
        mdata,
        **umap_params(filter_prevalence_df(mdata, "prop")),
        **filter_prevalence("prop")
    )


def ordinate_weighted(
    mdata: MuData,
    metadata_category: str,
    n_features=20
):
    ordinate(
        mdata,
        outputs_dest_key=f"umap_{metadata_category}",
        **umap_params(filter_diff_abund_df(mdata, "prop", metadata_category, n_features)),
        **filter_diff_abund("prop", metadata_category, n_features)
    )


def ordinate_features(mdata):

    ordinate(
        mdata,
        table=dict(
            data=dict(
                axis=1,
                tables=["prop.data"]
            )
        ),
        **umap_params(
            mdata.mod["prop"].to_df().T,
            metric="correlation"
        )
    )


def ordinate(
    mdata: MuData,
    outputs_dest_key="umap",
    metric="braycurtis",
    min_dist=0.2,
    n_components=2,
    n_neighbors=15,
    **kwargs
):
    """
    Run UMAP to ordinate the samples
    """
    # Run UMAP
    process.umap(
        mdata,
        umap_params=dict(
            metric=metric,
            min_dist=min_dist,
            n_components=n_components,
            n_neighbors=n_neighbors
        ),
        outputs_dest_key=outputs_dest_key,
        **kwargs
    )


def diff_abund(
    mdata: MuData,
    metadata_category,
    mod="prop",
):
    process.kruskal(
        mdata,
        outputs_dest_key=f"kruskal_{metadata_category}",
        table_grouping_axis=0,
        table_grouping_grouping_table=["Observation Metadata"],
        table_grouping_grouping_cname=metadata_category,
        table_grouping_grouping_label=metadata_category,
        **filter_prevalence(mod)
    )


def setup_config(df: DataFrame, config: Path):
    """Set up a config file based on the contents of the dataset."""

    # Pick the metadata column to use to compare samples
    cname, transform = pick_column(df)

    # Write the configuration file
    json.dump(
        [] if cname is None else [
            dict(
                dataset_name=(
                    config
                    .name
                    .replace(".config.json", "")
                    .replace("_", " ")
                ),
                cname=cname,
                label=cname.replace("_", " ").title() + (
                    f" ({transform})" if transform else ""
                ),
                n_features=20,
                transform=transform
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
        mdata = read_mdata(tsv)

        # Read the configuration file
        config_list = [
            Config(**cfg)
            for cfg in json.load(config.open())
        ]

        # Analyze each of the configured metadata categories
        for config_ix, config in enumerate(config_list):

            if config.cname is None:
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
