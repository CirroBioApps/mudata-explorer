from anndata import AnnData
from cirro import DataPortalDataset
from cirro.sdk.file import DataPortalFile
from muon import MuData
import pandas as pd
import streamlit as st
from typing import Dict, Tuple, Optional
from mudata_explorer.helpers.cirro_readers import util
from mudata_explorer.sdk import io, view

levels = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]


def read(
    dataset: DataPortalDataset
) -> Optional[MuData]:
    """Read datasets produced by curatedMetagenomicData."""

    # Get the parameters for formatting the MuData object
    mdata_params = _get_mdata_params(dataset)

    # If the parameters are not provided, return
    if mdata_params is None:
        return

    # Read the MuData object using those params
    with st.spinner("Reading in data"):
        mdata = _read_mdata(dataset, mdata_params)
        if mdata is None:
            return

    # Get the parameters for the analysis / display
    params = _get_params(mdata)
    if params is None:
        return

    if st.button("Load Dataset", key="load-from-cirro"):
        with st.spinner("Running Analysis"):
            _run_processes(mdata, params)
        with st.spinner("Configuring Displays"):
            _add_views(mdata, params)
        return mdata


def _get_mdata_params(dataset: DataPortalDataset) -> Optional[MuData]:
    # Get the list of files from the dataset
    files = util.list_files(dataset)

    # Filter to only those files which end with ".tsv" or ".csv"
    files = [
        f
        for f in files
        if f.name.endswith(".tsv") or f.name.endswith(".csv")
    ]

    # If there is more than one file, let the user select
    if len(files) > 1:
        selected_file = st.selectbox(
            "Select a file from the dataset",
            [f.name for f in files]
        )
        file = next(
            file
            for file in files
            if file.name == selected_file
        )
    else:
        file = files[0]

    # Select the taxonomic level to collapse by
    tax_level = st.selectbox(
        "Collapse by taxonomic level",
        levels,
        index=len(levels)-1
    )

    return dict(
        file_name=file.name,
        tax_level=tax_level
    )


def _read_mdata(
    dataset: DataPortalDataset,
    params: Dict[str, str]
) -> Optional[MuData]:

    # Get the list of files from the dataset
    files = util.list_files(dataset)

    # Get the specific file selected by the user
    file = next(
        f
        for f in files
        if f.name == params["file_name"]
    )

    # Read in the relative abundance table,
    # along with the taxonomic assignment for each feature
    # and the metadata for each sample
    obs, abund, var = _read_curated_metagenomic_data(file)

    # Collapse the abund and var by the indicated species
    abund, var = _collapse_by_taxon(abund, var, params["tax_level"])

    # Transform the abundances into proportions
    abund = abund.apply(lambda r: r / r.sum(), axis=1)

    # Make a MuData object
    mdata = io.build_mdata(
        dict(abund=AnnData(X=abund, var=var)),
        obs=obs
    )

    return mdata


def _read_curated_metagenomic_data(
    file: DataPortalFile
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Read the curatedMetagenomicData file format."""

    # Read in the table
    df: pd.DataFrame = file.read_csv(
        sep="\t" if file.name.endswith(".tsv") else ",",
        index_col=0
    )

    # Split it up into metadata and taxa
    obs = df.reindex(
        columns=[cname for cname in df.columns.values if '|' not in cname]
    )
    abund = df.reindex(
        columns=[cname for cname in df.columns.values if '|' in cname]
    )

    # Format the taxonomic levels for the variables
    levels_dict = {
        level.lower()[0]: level
        for level in levels
    }
    var = pd.DataFrame(
        [
            {
                levels_dict[org[0]]: org[3:]
                for org in taxon.split("|")
            }
            for taxon in abund.columns.values
        ],
        index=abund.columns.values
    )

    return obs, abund, var


def _collapse_by_taxon(
    abund: pd.DataFrame,
    var: pd.DataFrame,
    level: str
):

    assert level in var.columns.values, \
        f"Level '{level}' not found in taxonomic assignments."

    # Sum up the abundances for each taxon
    abund = (
        abund
        .T
        .groupby(var[level])
        .sum()
        .T
        .sort_index(axis=1)
        .rename(columns=_make_spaces)
    )

    # Collapse the taxonomic assignments
    var = (
        var
        .reindex(columns=levels[:levels.index(level)+1])
        .drop_duplicates()
        .set_index(level)
        .sort_index()
        .rename(index=_make_spaces)
        .apply(_make_spaces)
    )

    return abund, var


def _make_spaces(s: str):
    return s.replace("_", " ")


def _get_params(mdata: MuData) -> Optional[dict]:

    params = dict()

    # If there is sample metadata
    if mdata.obs.shape[1] > 0:

        # Let the user pick a column to compare samples by
        with st.container(border=1):
            st.write("""
**Compare Samples By:**

The user may optionally select a column from the sample metadata to
compare samples by.
The comparison will depend on whether the variable is continuous
(like age or weight) or categorical (discrete groups).
""")
            compare_by = st.selectbox(
                "Sample Metadata",
                ["< select a column >"] + list(mdata.obs.columns),
                index=0
            )

            if compare_by != "< select a column >":
                params["compare_by"] = compare_by

                # Get a readable label for this column
                st.write("**Variable Label**")
                params["label"] = st.text_input(
                    f"Label used in displays for '{compare_by}':",
                    compare_by.replace("_", " ")
                )

                params["is_categorical"] = util.ask_if_categorical(
                    params["compare_by"],
                    mdata.obs[params["compare_by"]],
                )

    with st.container(border=1):
        st.write("""**Top Features**

An extended set of figures will be displayed for the figures which are
the most strongly associated with the variable of interest.
The user may decide how many such variables to include in the display.
""")
        params["n_top_features"] = st.slider(
            "Number of features to display (with the lowest p-values)",
            min_value=1,
            max_value=100,
            value=10,
            step=1
        )

    with st.container(border=1):
        st.write("""**Sample Clustering**

The samples will be clustered into groups to identify
any subpopulations that share similar patterns of features.
                    """)
        params["leiden_res"] = st.slider(
            "Leiden Resolution - Higher values lead to more clusters",
            min_value=0.1,
            max_value=2.0,
            value=1.0,
            step=0.1
        )
    return params


def _run_processes(mdata: MuData, params: dict):

    # Summarize each feature
    with st.spinner("Summarizing organism abundances"):
        util.summarize_features(mdata, mod="abund")

    # Run UMAP
    with st.spinner("Running UMAP"):
        util.run_umap(
            mdata,
            mod="abund",
            metric="braycurtis",
            n_neighbors=15,
            dest_key="umap"
        )

    # Run leiden clustering on the samples
    with st.spinner("Running Leiden Clustering"):
        util.run_leiden(
            mdata,
            mod="abund",
            resolution=params["leiden_res"],
            metric="braycurtis",
            n_neighbors=15,
            dest_key="leiden"
        )

    # Run Kruskal to find the features which vary across clusters
    with st.spinner("Running Kruskal"):
        util.run_kruskal(
            mdata,
            table="abund.data",
            dest_key="kruskal_leiden",
            grouping_cname="leiden",
            grouping_table="Observation Metadata"
        )

    # If the user selected a metadata column
    # run Kruskal to find the features which vary across those groups
    if params.get("compare_by"):
        if params["is_categorical"]:
            with st.spinner("Running Kruskal"):
                util.run_kruskal(
                    mdata,
                    table="abund.data",
                    dest_key=f"kruskal_{params['compare_by']}",
                    grouping_cname=params["compare_by"],
                    grouping_table="Observation Metadata"
                )
        else:
            with st.spinner("Running Spearman"):
                util.run_spearman(
                    mdata,
                    table="abund.data",
                    dest_key=f"spearman_{params['compare_by']}",
                    comparitor_cname=params["compare_by"],
                    comparitor_table="Observation Metadata"
                )


def _add_views(mdata: MuData, params: dict):

    # Show the most abundant features
    util.add_boxplot_multi(
        mdata,
        mod="abund",
        title="Most Abundant Organisms",
        cols_query=dict(
            type="value",
            table=["abund.varm.summary_stats"],
            cname="mean_rank",
            expr="<=",
            value=str(params["n_top_features"])
        ),
        var_label="Species",
        val_label="Relative Abundance",
        legend="""The organisms which are most abundant in the dataset
are shown in the boxplot above.
The values displayed are relative which range from 0 to 1 -- calculated
as the proportion of reads from each sample which were assigned to
each organism."""
    )

    # Scatter plot with leiden clusters
    util.add_scatter(
        mdata,
        title="Map of Sample Similarity - Colored by Cluster",
        legend=f"""
Each point represents a single sample.
The samples are arranged in a two-dimensional space using the UMAP algorithm
such that samples with similar patterns of features are closer together.
The samples are colored based on unsupervised clustering with the
Leiden algorithm (resolution: {params['leiden_res']}).
""",
        table="abund.obsm.umap",
        axis=0,
        x="UMAP 1",
        xlabel="UMAP 1",
        y="UMAP 2",
        ylabel="UMAP 2",
        color_table="Observation Metadata",
        cname="leiden",
        label="Cluster (leiden)",
        is_categorical=True,
        scale="D3"
    )

    # Show the organisms that differentiate the leiden clusters
    view.plotly_category_summarize_values(
        mdata,
        **{
            "table": {
                "category": {
                    "category": {
                        "table": [
                            "Observation Metadata"
                        ],
                        "cname": "leiden",
                        "label": "Cluster (leiden)"
                    }
                },
                "data": {
                    "tables": [
                        "abund.data"
                    ],
                    "cols_query": {
                        "query": {
                            "cname": "rank",
                            "type": "value",
                            "table": [
                                "abund.varm.kruskal_leiden"
                            ],
                            "expr": "<=",
                            "value": str(params["n_top_features"])
                        }
                    },
                    "transforms": [
                        "zscores_cols"
                    ]
                }
            },
            "formatting": {
                "sort_by": "Values",
                "color": "None",
                "title": "Top Species by Cluster"
            }
        }
    )

    view.markdown(mdata, text="""The table above shows the top features which
most strongly vary between  the clusters identified by the Leiden algorithm.
The size of each point represents the mean abundance of the feature among the
samples which were assigned to a particular cluster.
To account for different scales of abundance, the values have been transformed
to z-scores.""")

    # Show the UMAP scatterplot, colored by the most differentiating feature
    top_feature = _get_top_feature(mdata, "kruskal_leiden")

    util.add_scatter(
        mdata,
        title=f"Map of Sample Similarity - Colored by {top_feature} Abundance",
        legend="""
The plot above shows how a single organism is distributed across samples.
Each point is colored by the relative abundance of the organism.
To create a similar display for any other organism, simply copy this
figure (select 'Edit Figures' on the left if the figure options are not
visible) and fill in the name of any organism to modify the display.""",
        table="abund.obsm.umap",
        axis=0,
        x="UMAP 1",
        xlabel="UMAP 1",
        y="UMAP 2",
        ylabel="UMAP 2",
        color_table="abund.data",
        cname=top_feature,
        label=top_feature,
        is_categorical=False,
        scale="bluered"
    )

    # If a metadata column was selected
    if params.get("compare_by"):

        # Show a UMAP, coloring by the variable of interest
        util.add_scatter(
            mdata,
            title=f"Map of Sample Similarity - Colored by {params['label']}",
            legend=f"""
Each point represents a single sample.
The samples are arranged in a two-dimensional space using the UMAP algorithm
such that samples with similar patterns of features are closer together.
The samples are colored based on the annotated value of '{params['label']}'.
    """,
            table="abund.obsm.umap",
            axis=0,
            x="UMAP 1",
            xlabel="UMAP 1",
            y="UMAP 2",
            ylabel="UMAP 2",
            color_table="Observation Metadata",
            cname=params["compare_by"],
            label=params["label"],
            is_categorical=params["is_categorical"],
            scale="D3" if params["is_categorical"] else "bluered",
        )

        # Comparison of the variable of interest across clusters
        title = f"Distribution of {params['label']} across Clusters"
        legend = f"""
    The distribution of the variable '{params['label']}'
    is shown across each of the clusters identified by the Leiden algorithm.
    """

        x_kwargs = dict(
            table="Observation Metadata",
            cname="leiden",
            label="Cluster (leiden)",
        )
        y_kwargs = dict(
            table="Observation Metadata",
            cname=params["compare_by"],
            label=params["label"]
        )

        if not params["is_categorical"]:

            # Display as a boxplot
            util.add_boxplot(
                mdata,
                title=title,
                legend=legend,
                **{
                    f"{prefix}_{kw}": val
                    for prefix, kwargs in [("x", x_kwargs), ("y", y_kwargs)]
                    for kw, val in kwargs.items()
                }
            )

        else:
            # For categorical variables, show a category
            # count which includes the leiden clusters
            util.add_category_count(
                mdata,
                title=title,
                legend=legend,
                **{
                    f"{prefix}_{kw}": val
                    for prefix, kwargs in [("x", x_kwargs), ("color", y_kwargs)]
                    for kw, val in kwargs.items()
                }
            )

        # Show the table of results
        stat_name = "kruskal" if params["is_categorical"] else "spearman"
        stats_table = f"abund.varm.{stat_name}_{params['compare_by']}"
        view.markdown(
            mdata,
            f"""**Comparison of Organisms by {params['label']}**

Each organism was individually compared to the variable '{params['label']}'.
Because the variable was marked as {'categorical' if params['is_categorical'] else 'continuous'},
the {stat_name.title()} test was used.
""")
        view.table(
            mdata,
            **{
                "options": {
                    "sort": {
                        "sort_by": {
                            "table": [stats_table],
                            "cname": "pvalue",
                            "label": "pvalue"
                        },
                        "axis": 1
                    }
                },
                "data": {
                    "table": {
                        "axis": 1,
                        "tables": [stats_table]
                    }
                }
            }
        )

        # Get the top organism
        top_feature = _get_top_feature(mdata, f"{stat_name}_{params['compare_by']}")

        # Now show the specific organisms which vary by metadata
        # If the variable is categorical
        if params["is_categorical"]:

            # Dotplot with organisms associated with the variable of interest
            view.plotly_category_summarize_values(
                mdata,
                **{
                    "table": {
                        "category": {
                            "category": {
                                "table": [
                                    "Observation Metadata"
                                ],
                                "cname": params["compare_by"],
                                "label": params["label"]
                            }
                        },
                        "data": {
                            "tables": [
                                "abund.data"
                            ],
                            "cols_query": {
                                "query": {
                                    "cname": "rank",
                                    "type": "value",
                                    "table": [
                                        f"abund.varm.kruskal_{params['compare_by']}"
                                    ],
                                    "expr": "<=",
                                    "value": str(params["n_top_features"])
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

            # Show a multi boxplot with the top organisms
            view.plotly_box_multiple(
                mdata,
                **{
                    "table": {
                        "category": {
                            "category": {
                                "table": [
                                    "Observation Metadata"
                                ],
                                "cname": "ici_response",
                                "label": "ici_response"
                            },
                            "rows_query": {
                                "query": {
                                    "cname": "",
                                    "type": "value"
                                }
                            }
                        },
                        "data": {
                            "tables": [
                                "abund.data"
                            ],
                            "cols_query": {
                                "query": {
                                    "cname": "rank",
                                    "type": "value",
                                    "table": [
                                        "abund.varm.kruskal_ici_response"
                                    ],
                                    "expr": "<=",
                                    "value": "10"
                                }
                            },
                            "rows_query": {
                                "query": {
                                    "cname": "",
                                    "type": "value"
                                }
                            },
                            "transforms": []
                        }
                    },
                    "category_options": {
                        "axis": "Color"
                    },
                    "display_options": {
                        "outliers": {
                            "enabled": True
                        }
                    },
                    "variable_options": {
                        "log_values": False,
                        "axis": "X-Axis"
                    }
                }
            )

        # If the variable is continuous
        else:

            # Show a scatterplot with the top organism
            util.add_scatter(
                mdata,
                title=f"Abundance of {top_feature} by {params['label']}",
                legend=f"""As an example of the relationship between a single
organism and the {params['label']} variable, the scatterplot above shows the
abundance of '{top_feature}' across the samples.
The color of each point indicates the cluster assigned to the sample
by the Leiden algorithm.
""",
                table=None,
                xtable="Observation Metadata",
                x=params["compare_by"],
                xlabel=params["label"],
                ytable="abund.data",
                y=top_feature,
                ylabel=top_feature,
                axis=0,
                color_table="Observation Metadata",
                cname="leiden",
                label="Cluster",
                is_categorical=True,
                scale="D3"
            )


def _get_top_feature(mdata: MuData, table: str) -> str:
    return (
        mdata
        .mod["abund"]
        .varm[table]
        .sort_values("rank")
        .index[0]
    )
