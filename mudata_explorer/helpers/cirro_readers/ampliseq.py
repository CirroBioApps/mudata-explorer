from anndata import AnnData
from cirro import DataPortalDataset
from cirro.sdk.file import DataPortalFile
from muon import MuData
import pandas as pd
import streamlit as st
from typing import Tuple, List, Optional
from mudata_explorer.helpers.cirro_readers import util
from mudata_explorer.sdk import io, view

levels = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Species_exact"]


def read(
    dataset: DataPortalDataset
) -> Optional[MuData]:
    """Read datasets produced by nf-core/ampliseq."""

    mdata = _read_mdata(dataset)
    if mdata is None:
        return

    params = _get_params(mdata)
    if params is None:
        return

    if st.button("Load Dataset", key="load-from-cirro"):
        _run_processes(mdata, params)
        with st.spinner("Configuring Displays"):
            _add_views(mdata, params)
        return mdata


def _read_mdata(dataset: DataPortalDataset) -> Optional[MuData]:
    # Get the list of files from the dataset
    files = util.list_files(dataset)

    # Read in the relative abundance table,
    # along with the taxonomic assignment for each feature
    with st.container(border=1):
        tax, abund = _read_rel_abund(files)

    # Get the sample metadata
    with st.container(border=1):
        sample_meta = _read_sample_meta(dataset, abund)

    if sample_meta is None:
        return

    # For this dataset, any dashes in sample names
    # are replaced with underscores
    sample_meta = sample_meta.rename(
        index=lambda i: i.replace("-", "_")
    )

    # Make a MuData object
    mdata = io.build_mdata(
        dict(
            abund=AnnData(
                X=abund.sort_index().T,
                var=tax.sort_index()
            )
        ),
        obs=sample_meta.reindex(abund.columns)
    )

    return mdata


def _read_rel_abund(
    files: List[DataPortalFile]
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Read the relative abundance table from a dataset."""

    # Filter to only the relative abundance tables
    abund_prefix = "data/qiime2/rel_abundance_tables/rel-table-ASV_with-"
    abund_suffix = "-tax.tsv"
    abund_files = [
        f
        for f in files
        if f.name.startswith(abund_prefix) and f.name.endswith(abund_suffix)
    ]

    # If there are no relative abundance tables
    if len(abund_files) == 0:
        st.error("No relative abundance tables found.")
        return
    # If there is more than one relative abundance table
    elif len(abund_files) > 1:
        # Pick the file from the dataset
        selected_abund = st.selectbox(
            "Use taxonomy assigned by:",
            [f.name[len(abund_prefix):-len(abund_suffix)] for f in abund_files]
        )
        abund_file = next(
            file
            for file in abund_files
            if file.name == f"{abund_prefix}{selected_abund}{abund_suffix}"
        )
    else:
        abund_file = abund_files[0]

    # Read in the table
    df = abund_file.read_csv(sep="\t", index_col=0)

    # Drop the unused columns
    df = (
        df
        .rename(columns=dict(confidence="Confidence"))
        .drop(columns=["sequence", "Confidence"] + [
            col for col in ["Domain"]
            if col in df.columns.values
        ])
    )

    # Parse the taxonomy
    df = _parse_taxonomy(df)

    # Fix the species names
    df = _fix_species(df)

    # Collapse by taxon
    return _collapse_by_taxon(df)


def _collapse_by_taxon(df: pd.DataFrame):
    filt_levels = [level for level in levels if level in df.columns]
    assert len(filt_levels) > 0

    # Let the user select which taxon to collapse by
    level = st.selectbox(
        "Label organisms to the level of:",
        filt_levels,
        index=max(0, len(filt_levels)-2)
    )

    # If the ASV isn't classified at that level, fill in with the
    # next highest level
    df = df.assign(
        _group_level=df.apply(
            lambda r: _pick_level(r, level),
            axis=1
        )
    )

    # Get the taxonomic assignment for each taxon
    tax = (
        df.reindex(
            columns=filt_levels[:filt_levels.index(level)] + ["_group_level"]
        )
        .groupby("_group_level")
        .head(1)
        .set_index("_group_level")
        .rename_axis(index=level)
    )

    abund = (
        df
        .drop(columns=filt_levels)
        .groupby("_group_level")
        .sum()
        .rename_axis(index=level)
    )

    return tax, abund


def _pick_level(row: pd.Series, level: str):

    return (
        row[level]
        if not pd.isnull(row.get(level))
        else next(
            f"{row[fallback_level]} ({fallback_level})"
            for fallback_level in levels[::-1]
            if not pd.isnull(row.get(fallback_level))
        )
    )


@st.cache_data
def _parse_taxonomy(df: pd.DataFrame):
    if "Taxon" not in df.columns:
        return df

    return pd.concat(
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
    }


def _fix_species(df: pd.DataFrame):
    for kw in ["Species", "Species_exact"]:
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


def _read_sample_meta(
    dataset: DataPortalDataset,
    abund: pd.DataFrame
) -> pd.DataFrame:

    # Read the samplesheet saved with the dataset
    meta = util.sample_metadata(dataset)

    # For this process, any samples with numeric names have "Sample_" prepended
    meta = meta.rename(
        index=lambda i: f"Sample_{i}" if i[0] in map(str, range(10)) else i
    )

    st.write("""
**Sample Metadata**

The table below shows the metadata associated with each sample in the dataset.
To modify or augment this information, download as a CSV and upload a new copy.
""")

    # Let the user view the metadata
    st.write(meta.reset_index())

    # Provide a button to download as CSV
    st.download_button(
        label="Download sample metadata",
        data=meta.to_csv(),
        file_name="sample_metadata.csv",
        mime="text/csv"
    )

    # Let the user upload a new samplesheet
    uploaded_file = st.file_uploader("Upload a new spreadsheet (in CSV or XLSX format)")

    if uploaded_file:
        if uploaded_file.name.endswith(".csv"):
            new_meta = pd.read_csv(uploaded_file, index_col=0)
        elif uploaded_file.name.endswith(".xlsx"):
            new_meta = pd.read_excel(uploaded_file, index_col=0)
        else:
            st.error("Please upload a CSV or XLSX file.")
            return

        # Subset to only those rows which are shared
        n_intersection = new_meta.index.intersection(meta.index).shape[0]
        st.write(f"The uploaded file has {n_intersection:,} valid sample labels.")
        if n_intersection > 0:

            # Drop any columns which are entirely missing
            new_meta = new_meta.drop(
                columns=[
                    cname
                    for cname in new_meta.columns
                    if new_meta[cname].isnull().all()
                ]
            )

            # Show the user the updated metadata
            new_meta = new_meta.reindex(meta.index)
            st.dataframe(new_meta)
            return new_meta

    else:
        return meta


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


def _add_views(mdata: MuData, params: dict):

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
    view.markdown(mdata, text="**Top Features by Cluster**")

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
                "color": "None"
            }
        }
    )

    view.markdown(mdata, text="""The table above shows the top features which
are most strongly associated with the clusters identified by the Leiden
algorithm.
The size of each point represents the mean abundance of the feature among the
samples which were assigned to a particular cluster.
To account for different scales of abundance, the values have been transformed
to z-scores.""")

    # If a metadata column was selected
    if params.get("compare_by"):

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
