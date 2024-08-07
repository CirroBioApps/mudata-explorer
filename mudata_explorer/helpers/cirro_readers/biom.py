from tempfile import TemporaryDirectory
from biom import load_table
from anndata import AnnData
from cirro import DataPortalDataset
from cirro.sdk.file import DataPortalFile
from muon import MuData
import pandas as pd
import streamlit as st
from typing import List, Optional
from mudata_explorer.helpers.cirro_readers import util
from mudata_explorer.sdk import io, view


def read(
    dataset: DataPortalDataset
) -> Optional[MuData]:
    """Read BIOM datasets."""

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

    # Read in a BIOM file as AnnData
    adata = _read_biom_as_anndata(files)
    if adata is None:
        return

    # Pick a taxonomic level to collapse by
    with st.container(border=1):
        adata = _collapse_by_taxon(adata)

    # Get the sample metadata
    with st.container(border=1):
        adata = _read_sample_meta(adata)

    # Make a MuData object
    mdata = io.build_mdata(
        dict(
            abund=AnnData(
                X=adata.to_df(),
                var=adata.var
            )
        ),
        obs=adata.obs
    )

    return mdata


def _read_biom_as_anndata(
    files: List[DataPortalFile]
) -> Optional[AnnData]:
    """Read the BIOM file from a dataset."""

    biom_file = util.find_file_by_extension(
        files,
        prefix="data/",
        suffix=".biom",
        selectbox_label="Use BIOM file:",
        none_found_msg="No BIOM files found."
    )
    if biom_file is None:
        return

    # Download the file and read it
    with st.spinner("Reading BIOM file"):
        with TemporaryDirectory() as tmpdir:
            biom_file.download(tmpdir)
            table = load_table(f"{tmpdir}/{biom_file.name}")
    return table.to_anndata()


def _collapse_by_taxon(adata: AnnData) -> AnnData:
    # The taxonomy is in the var table. However, it is not labelled.
    filt_levels = list(adata.var.columns.values)
    assert len(filt_levels) > 0

    # Let the user select which taxon to collapse by
    level = st.selectbox(
        "Label organisms to the level of:",
        filt_levels,
        index=max(0, len(filt_levels)-2)
    )

    # Show the user the most common taxa at this level
    st.dataframe(
        adata
        .var[level]
        .fillna("Unclassified")
        .value_counts()
        .reset_index()
        .rename(columns=dict(count="Number of ASVs")),
        hide_index=True
    )

    # If the ASV isn't classified at that level, fill in with the
    # next highest level
    adata.var["_group_level"] = adata.var.apply(
        lambda r: _pick_level(r, filt_levels.index(level)),
        axis=1
    )

    # Sum up the abundances for each taxon
    abund = (
        adata
        .to_df()
        .T
        .groupby(adata.var["_group_level"])
        .sum()
        .sort_index()
        .rename_axis(index=level)
        .T
    )

    # Get the taxonomic assignment for each taxon
    tax = (
        adata.var
        .rename(
            columns={
                f"taxonomy_{ix}": label
                for ix, label in enumerate(filt_levels)
            }
        )
        .reindex(columns=filt_levels[:filt_levels.index(level)] + ["_group_level"])
        .groupby("_group_level")
        .head(1)
        .set_index("_group_level")
        .sort_index()
        .rename_axis(index=level)
        .reindex(index=abund.columns.values)
    )

    # Return an AnnData object
    return AnnData(X=abund, var=tax, obs=adata.obs)


def _pick_level(row: pd.Series, pick_ix: int):

    # Walk from each of the levels back to the first
    for ix in range(pick_ix, -1, -1):

        # If there is a value at this level
        if not pd.isnull(row.iloc[ix]):
            # If we're at a higher level than intended,
            # annotate the rank
            if ix < pick_ix:
                return f"{row.iloc[ix]} ({levels[ix]})"
            # If we are at the desired level, no need to label
            else:
                return row.iloc[ix]

    # If none of the levels have values
    return "Unclassified"


def _read_sample_meta(adata: AnnData) -> AnnData:

    st.write("""
**Sample Metadata**

The table below shows the metadata associated with each sample in the dataset.
To modify or augment this information, download as a CSV and upload a new copy.
""")

    # Let the user view the metadata
    st.write(adata.obs.reset_index())

    # Provide a button to download as CSV
    st.download_button(
        label="Download sample metadata",
        data=adata.obs.to_csv(),
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
        n_intersection = new_meta.index.intersection(adata.obs.index).shape[0]
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
            new_meta = new_meta.reindex(adata.obs.index)
            st.dataframe(new_meta)
            adata.obs = new_meta

    return adata


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
