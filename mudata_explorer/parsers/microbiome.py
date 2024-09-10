from anndata import AnnData
from copy import copy
from dataclasses import dataclass
from mudata_explorer.sdk import io, view
from mudata_explorer.helpers.cirro_readers import util
from mudata_explorer.sdk import view, io
from muon import MuData
import pandas as pd
from typing import Optional
import streamlit as st


@dataclass
class MicrobiomeParams:
    """Configuration for the analysis."""
    # Name of the dataset, to be used in the title
    dataset_name: str = "Microbiome Analysis"
    # Column name to use for comparison
    compare_by: Optional[str] = None
    # Label to use when describing those comparison groups
    label: Optional[str] = None
    # Number of features to use in the analysis and plots
    n_top_features: int = 20
    # Whether the comparison variable is categorical
    is_categorical: Optional[str] = None
    # Resolution for the Leiden clustering
    leiden_res: float = 1.0


def _top_features(
    mdata: MuData,
    n_top_features: int,
    varm_key="summary_stats",
    cname="mean",
    ascending=False
):
    return list(
        mdata
        .mod["abund"]
        .varm[varm_key]
        .sort_values(
            by=cname,
            ascending=ascending
        )
        .head(n_top_features)
        .index.values
    )


def _top_features_filter_cols(
    mdata: MuData,
    n_top_features: int,
    varm_key="summary_stats",
    cname="mean",
    ascending=False
): 
    features = _top_features(
        mdata,
        n_top_features,
        varm_key=varm_key,
        cname=cname,
        ascending=ascending
    )

    return dict(
        tables_value=["abund.data"],
        type_value="index",
        expr_value="in",
        value_enum_value=features,
        value_enum_sidebar=True
    )


def parse_adata(adata: AnnData, groupby_var=False) -> Optional[MuData]:

    mdata = _parse_adata(adata, groupby_var=groupby_var)

    params = _get_params(mdata)
    if params is None:
        return

    if st.button("Load Dataset", key="load-from-cirro"):
        _run_processes(mdata, params)
        with st.spinner("Configuring Displays"):
            _add_views(mdata, params)
        return mdata


def _parse_adata(adata: AnnData, groupby_var=False) -> MuData:

    # Filter out any samples with 0 counts
    adata = adata[adata.to_df().sum(axis=1) > 0]
    assert adata.shape[0] > 0, "No samples with non-zero counts."

    # Make sure that every row sums to 1
    adata.X = adata.to_df().apply(lambda r: r / r.sum(), axis=1).values

    if groupby_var:
        # Pick a taxonomic level to collapse by
        with st.container(border=1):
            adata = _collapse_by_taxon(adata)

    # Get the sample metadata
    with st.container(border=1):
        adata = _read_sample_meta(adata)

    # Optionally filter samples by sample metadata
    with st.container(border=1):
        adata = _filter_samples(adata)

    # Coerce objects to str in the metadata table
    adata.obs = adata.obs.astype({
        cname: "str"
        for cname, dtype in adata.obs.dtypes.items()
        if dtype == "object"
    })

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

    # Show the user the most abundant taxa at this level
    st.dataframe(
        abund
        .mean()
        .sort_values(ascending=False)
        .reset_index()
        .rename(columns={0: "Average Relative Abundance"}),
        hide_index=True
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


def _filter_samples(adata: AnnData):

    st.markdown("""**Filter Samples**

Optionally filter the set of samples which are included in the analysis
""")
    
    # Let the user select the approach for filtering samples
    filter_method = st.selectbox(
        "Filter Samples",
        [
            "Query String",
            "Select Samples to Include",
            "Select Samples to Exclude"
        ],
        index=None,
        placeholder="Select a method"
    )
    if filter_method == "Query String":
        return _filter_samples_query(adata)
    elif filter_method == "Select Samples to Include":
        return _filter_samples_select(adata, include=True)
    elif filter_method == "Select Samples to Exclude":
        return _filter_samples_select(adata, include=False)
    else:
        return adata


def _filter_samples_query(adata: AnnData) -> AnnData:
    ntot = adata.shape[0]
    
    query_string = st.text_input(
        "Query String",
        placeholder="e.g. group == 'A' or timepoint == 1"
    )
    if query_string is not None and len(query_string) > 0:
        try:
            filtered = adata[
                adata.obs.query(query_string).index.values
            ]
        except Exception as e:
            st.exception(e)

        if filtered.shape[0] > 0:

            st.markdown(f"Including {filtered.shape[0]:,} / {ntot:,} samples")
            return filtered
        
        else:
            st.write("No samples matching the query")

    else:
        st.markdown(f"Including all {ntot:,} samples")

    return adata


def _filter_samples_select(adata: AnnData, include: bool) -> AnnData:
    ntot = adata.shape[0]

    # Let the user select samples to include or exclude
    selected = st.multiselect(
        f"Select Samples to {'Include' if include else 'Exclude'}",
        adata.obs.index.values,
        adata.obs.index.values if include else []
    )

    if len(selected) > 0 and len(selected) < ntot:

        if include:
            filtered = adata[selected]
            st.markdown(f"Including {filtered.shape[0]:,} / {ntot:,} samples")
        else:
            filtered = adata[~adata.obs.index.isin(selected)]
            st.markdown(f"Including {filtered.shape[0]:,} / {ntot:,} samples")

        return filtered

    else:
        st.markdown(f"Including all {ntot:,} samples")
        return adata


def _get_params(mdata: MuData) -> Optional[MicrobiomeParams]:

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

    params["dataset_name"] = st.text_input(
        "Dataset Name",
        "Microbiome Dataset"
    )

    return MicrobiomeParams(**params)


def _run_processes(mdata: MuData, params: MicrobiomeParams):

    # Get summary metrics for each feature
    with st.spinner("Calculating Feature Metrics"):
        util.summarize_features(mdata, mod="abund")

    # Calculate alpha diversity for each sample
    with st.spinner("Calculating Alpha Diversity"):
        util.shannon_diversity(mdata, mod="abund")

    # Run PCA
    with st.spinner("Running PCA"):
        util.run_pca(
            mdata,
            mod="abund",
            dest_key="pca"
        )

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
            resolution=params.leiden_res,
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
    if params.compare_by:

        # If the metadata is categorical
        if params.is_categorical:
                
            # Run Kruskal to find the features which vary across groups
            with st.spinner("Running Kruskal"):
                util.run_kruskal(
                    mdata,
                    table="abund.data",
                    dest_key=f"kruskal_{params.compare_by}",
                    grouping_cname=params.compare_by,
                    grouping_table="Observation Metadata"
                )

        # If the metadata is continuous
        else:

            # Run Spearman to find the
            # features which correlate with the variable
            with st.spinner("Running Spearman"):
                util.run_spearman(
                    mdata,
                    table="abund.data",
                    dest_key=f"spearman_{params.compare_by}",
                    comparitor_cname=params.compare_by,
                    comparitor_table="Observation Metadata"
                )


def _view_stacked_bars(mdata: MuData, params: MicrobiomeParams):

    # Show some stacked bars
    util.add_stacked_bars(
        mdata,
        title="Taxonomic Composition",
        yaxis_title="Proportion of Reads",
        table="abund.data",
        # Only show the top features
        features=_top_features(mdata, params.n_top_features),
        category_cname=params.compare_by if params.is_categorical else None,
        category_label=params.label,
        legend="""The taxonomic composition of this group of samples is shown using
a stacked bar plot. Each bar represents a single sample, and the height of the
bar is proportional to the relative abundance of each taxon in that sample.
The relative abundance of all organisms sums to 1 for each sample.
""" + (f"The samples are colored based on the annotated value of '{params.label}'." if params.is_categorical else "")
    )


def _view_shannon_diversity(mdata: MuData, params: MicrobiomeParams):

    # Set up the basic kwargs for the histogram
    kwargs = {
        "data": {
            "sidebar": False,
            "axis": {"value": 0, "sidebar": False},
            "columns": {
                "value": {
                    "sidebar": False,
                    "table": {
                        "value": "Observation Metadata",
                        "sidebar": False
                    },
                    "cname": {
                        "value": "shannon",
                        "sidebar": False
                    },
                    "label": {
                        "value": "Shannon Diversity Index",
                        "sidebar": False
                    },
                    "colorscale": False
                },
                "grouping": {
                    "enabled": {
                        "value": False,
                        "sidebar": False
                    },
                    "sidebar": False
                }
            }
        },
        "scale_options": {
            "log_y": {
                "value": None,
                "sidebar": False
            }
        },
        "formatting": {
            "title": {
                "value": "Distribution of Alpha Diversity (Shannon Index)",
                "sidebar": True
            },
            "legend": {
                "value": """The Shannon Diversity Index is a measure of the diversity
of organisms in a sample. It takes into account both the number of organisms
and the evenness of their distribution.
The height of each bar in this frequency histogram depicts the number of samples
which have a particular value of the Shannon Diversity Index (shown on the x-axis).
""",
                "sidebar": True
            },
            "nbins": {
                "value": 30,
                "sidebar": True
            }
        },
        "statistics": {
            "compare_groups": {
                "value": "Disabled",
                "sidebar": False
            }
        }
    }

    # Show the histogram
    view.plotly_histogram(
        mdata,
        **copy(kwargs)
    )

    # If the user selected a categorical metadata column to compare between
    if params.compare_by and params.is_categorical:

        # Enable the grouping
        kwargs["data"]["columns"]["grouping"] = {
            "enabled": {
                "value": True,
                "sidebar": False
            },
            "sidebar": False,
            "table": {
                "value": "Observation Metadata",
                "sidebar": False
            },
            "cname": {
                "value": params.compare_by,
                "sidebar": False
            },
            "label": {
                "value": params.label,
                "sidebar": False
            },
            "scale": {
                "value": "D3",
                "sidebar": False
            },
            "colorscale": True,
            "is_categorical": {
                "value": True,
                "sidebar": False
            }
        }
        # Update the title
        kwargs["formatting"]["title"]["value"] = f"Distribution of Alpha Diversity by {params.label}"
        # Update the legend
        kwargs["formatting"]["legend"]["value"] += f"The samples are colored based on the annotated value of '{params.label}'."
        # Compute stats between the groups
        kwargs["statistics"]["compare_groups"]["value"] = "Kruskal-Wallis"

        # Add that histogram as well
        view.plotly_histogram(
            mdata,
            **copy(kwargs)
        )


def _view_most_abundant_boxplot(mdata: MuData, params: MicrobiomeParams):

    # Show the most abundant organisms
    view.plotly_box_multiple(
        mdata,
        table_category_enabled_value=False,
        table_data_tables_value=["abund.data"],
        table_data_filter_cols=_top_features_filter_cols(mdata, params.n_top_features),
        variable_options_axis_value="X-Axis",
        variable_options_log_values_value=False,
        display_options_title_value="High Abundance Organisms",
        display_options_var_label_value="Organisms",
        display_options_var_label_sidebar=True,
        display_options_val_label_value="Relative Abundance",
        display_options_val_label_sidebar=True,
        display_options_legend_value="""The box plot above shows the relative abundance of the most
abundant organisms in the dataset. Each box represents the distribution of
abundances across the samples for a single organism."""
    )

    # If the user selected a categorical metadata column to compare between
    if params.compare_by and params.is_categorical:
            
        # Show the most abundant organisms by group
        view.plotly_box_multiple(
            mdata,
            table_category_enabled_value=True,
            table_category_columns_category_table_value="Observation Metadata",
            table_category_columns_category_cname_value=params.compare_by,
            table_category_columns_category_label_value=params.label,
            table_data_tables_value=["abund.data"],
            table_data_filter_cols=_top_features_filter_cols(mdata, params.n_top_features),
            variable_options_axis_value="X-Axis",
            variable_options_log_values_value=False,
            display_options_title_value=f"High Abundance Organisms by {params.label}",
            display_options_var_label_value="Organisms",
            display_options_var_label_sidebar=True,
            display_options_val_label_value="Relative Abundance",
            display_options_val_label_sidebar=True,
            category_options_axis_value="Color",
            category_options_sort_by_value="Mean",
            display_options_legend_value=f"""The box plot above shows the relative abundance of the most
abundant organisms in the dataset. Each box represents the distribution of
abundances across the samples for a single organism. The boxes are colored
based on the annotated value of '{params.label}'."""
        )


def _view_pca(mdata: MuData, params: MicrobiomeParams):
    pc1_cname = mdata.mod["abund"].obsm["pca"].columns[0]
    pc2_cname = mdata.mod["abund"].obsm["pca"].columns[1]

    # Get the organisms with the largest loadings on PC1 and PC2
    top_orgs = list(
        mdata.mod["abund"]
        .varm["pca"]
        .reindex(columns=[pc1_cname, pc2_cname])
        .abs()
        .max(axis=1)
        .sort_values(ascending=False)
        .head(params.n_top_features)
        .index.values
    )

    # Show the PCA, coloring the points if a metadata column was selected
    view.plotly_scatter_vectors(
        mdata,
        data_axis_value=0,
        data_columns_x_table_value="abund.obsm.pca",
        data_columns_y_table_value="abund.obsm.pca",
        data_columns_x_cname_value=pc1_cname,
        data_columns_x_label_value=pc1_cname,
        data_columns_y_cname_value=pc2_cname,
        data_columns_y_label_value=pc1_cname,
        data_columns_size_enabled_value=False,
        data_columns_color_enabled_value=params.compare_by is not None,
        data_columns_color_table_value="Observation Metadata",
        data_columns_color_cname_value=params.compare_by,
        data_columns_color_label_value=params.label,
        data_columns_color_colorscale=("D3" if params.is_categorical else "bluered"),
        data_columns_color_is_categorical_value=params.is_categorical,
        vectors_axis_value=1,
        vectors_columns_x_table_value="abund.varm.pca",
        vectors_columns_y_table_value="abund.varm.pca",
        vectors_columns_x_cname_value=pc1_cname,
        vectors_columns_x_label_value=pc1_cname,
        vectors_columns_y_cname_value=pc2_cname,
        vectors_columns_y_label_value=pc1_cname,
        vectors_columns_label_enabled_value=False,
        vectors_filter_rows_type_value="index",
        vectors_filter_rows_tables_value=["abund.varm.pca"],
        vectors_filter_rows_expr_value="in",
        vectors_filter_rows_value_enum_value=top_orgs,
        vectors_filter_rows_value_enum_sidebar=True,
        formatting_title_value="Beta Diversity PCA",
        formatting_legend_value=f"""Each point represents a single sample.
The samples are arranged in a two-dimensional space using the PCA algorithm
such that samples with similar patterns of features are closer together.
The color of each point represents the annotated value of '{params.label}'."""
    )


def _view_umap(mdata: MuData, params: MicrobiomeParams):
    # If a metadata column was selected
    if params.compare_by:

        # Show a UMAP, coloring by the variable of interest
        util.add_scatter(
            mdata,
            title=f"Map of Sample Similarity - Colored by {params.label}",
            legend=f"""
Each point represents a single sample.
The samples are arranged in a two-dimensional space using the UMAP algorithm
such that samples with similar patterns of features are closer together.
The samples are colored based on the annotated value of '{params.label}'.
    """,
            table="abund.obsm.umap",
            axis=0,
            x="UMAP 1",
            xlabel="UMAP 1",
            y="UMAP 2",
            ylabel="UMAP 2",
            color_table="Observation Metadata",
            cname=params.compare_by,
            label=params.label,
            is_categorical=params.is_categorical,
            scale="D3" if params.is_categorical else "bluered",
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
Leiden algorithm (resolution: {params.leiden_res}).
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


def _view_orgs_by_cluster(mdata: MuData, params: MicrobiomeParams):

    # Show the most variable organisms by cluster
    view.plotly_category_summarize_values(
        mdata,
        table_category_columns_category_table_value="Observation Metadata",
        table_category_columns_category_cname_value="leiden",
        table_category_columns_category_label_value="Cluster (leiden)",
        table_data_tables_value=["abund.data"],
        table_data_filter_cols=_top_features_filter_cols(
            mdata,
            params.n_top_features,
            varm_key="summary_stats",
            cname="std",
            ascending=False
        ),
        table_data_transforms_value=["zscores_cols"],
        formatting_sort_by_value="Values",
        formatting_color_value="None",
        formatting_title_value="Top Organisms by Cluster",
        formatting_legend_value="""The table above shows the top features which
have the highest standard deviation of abundances across the samples.
The abundance of each organism is shown for each of the clusters identified by the
Leiden algorithm.
The size of each point represents the mean abundance of the feature among the
samples which were assigned to a particular cluster.
To account for different scales of abundance, the values have been transformed
to z-scores.""",
        formatting_legend_sidebar=True
    )


def _view_grouping_across_clusters(mdata: MuData, params: MicrobiomeParams):

    # Comparison of the variable of interest across clusters
    title = f"Distribution of {params.label} across Clusters"
    legend = f"""The distribution of the variable '{params.label}'
is shown across each of the clusters identified by the Leiden algorithm.
"""

    x_kwargs = dict(
        table="Observation Metadata",
        cname="leiden",
        label="Cluster (leiden)",
    )
    y_kwargs = dict(
        table="Observation Metadata",
        cname=params.compare_by,
        label=params.label
    )

    if not params.is_categorical:

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


def _setup_stats_table_name(params: MicrobiomeParams):
    if params.is_categorical:
        stats_table = f"abund.varm.kruskal_{params.compare_by}"
        stats_name = "Kruskal-Wallis H-test"
    else:
        stats_table = f"abund.varm.spearman_{params.compare_by}"
        stats_name = "Spearman Rank Correlation"
    return stats_table, stats_name


def _view_orgs_assoc_with_clusters(mdata: MuData, params: MicrobiomeParams):
    """Summarize the association of individual organisms with the variable"""

    stats_table, stats_name = _setup_stats_table_name(params)

    view.table(
        mdata,
        options_sort_axis_value=1,
        options_sort_columns_sort_by_table_value=[stats_table],
        options_sort_columns_sort_by_cname_value="pvalue",
        options_sort_columns_sort_by_label_value="pvalue",
        data_table_axis_value=1,
        data_table_tables_value=[stats_table],
        data_table_filter_cols_tables_value=[stats_table],
        data_table_filter_cols_type_value="index",
        data_table_filter_cols_expr_value="in",
        data_table_filter_cols_cname_value="",
        data_table_filter_cols_value_enum_value=[
            "statistic",
            "pvalue",
            "neg_log10_pvalue",
            "mean"
        ],
        data_table_filter_cols_value_enum_sidebar=True,
        display_options_title_value=f"Top Organisms Associated with {params.label}",
        display_options_legend_value=f"""
The table above shows the top features which are most strongly associated with
the variable of interest ({params.label}).
{stats_name} was used to identify the features which are most strongly
associated with the variable."""
    )

    view.plotly_scatter(
        mdata,
        data_axis_value=1,
        data_columns_color_enabled_value=False,
        data_columns_size_enabled_value=False,
        data_columns_x_table_value=stats_table,
        data_columns_x_cname_value="mean",
        data_columns_x_label_value="Mean Abundance",
        data_columns_x_label_sidebar=True,
        data_columns_y_table_value=stats_table,
        data_columns_y_cname_value="neg_log10_pvalue",
        data_columns_y_label_value="-log10(p-value)",
        data_columns_y_label_sidebar=True,
        formatting_title_value=f"Association with {params.label} ({stats_name})",
        formatting_title_sidebar=True,
        scale_options_log_x_value=True,
        scale_options_log_x_sidebar=True,
        scale_options_log_y_sidebar=False,
        formatting_legend_value=f"""Each point represents a single microbe.
The vertical position (y-axis) of each point shows the degree of association using the -log10(p-value).
The x-axis shows the mean abundance of the microbe among the samples."""
    )


def _view_top_assoc_org(mdata: MuData, params: MicrobiomeParams):

    stats_table, _ = _setup_stats_table_name(params)

    # Plot the most strongly associated organism
    org = mdata.mod["abund"].varm[stats_table.split(".")[-1]]["neg_log10_pvalue"].idxmax()
    if params.is_categorical:
        util.add_boxplot(
            mdata,
            title="",
            legend=f"""The distribution of abundances for a single organism are shown
across the different groups identified by '{params.label}'.""",
            x_table="Observation Metadata",
            x_cname=params.compare_by,
            x_label=params.label,
            y_table="abund.data",
            y_cname=org,
            y_label=org
        )
    else:
        util.add_scatter(
            mdata,
            title="",
            legend=f"""The abundance of a single organism is shown across the samples.
The x-axis shows the value of {params.label} for each sample.
The y-axis shows the relative abundance of a single organism.""",
            x=params.compare_by,
            xlabel=params.label,
            y=org,
            ylabel=org,
            table=None,
            axis=0,
            color_table=None,
            cname=None,
            label=None,
            is_categorical=False,
            scale=None,
            xtable="Observation Metadata",
            ytable="abund.data"
        )


def _add_views(mdata: MuData, params: MicrobiomeParams):

    if params.dataset_name:
        view.markdown(mdata, f"### {params.dataset_name}", text_sidebar=True)

    _view_stacked_bars(mdata, params)

    _view_shannon_diversity(mdata, params)

    _view_most_abundant_boxplot(mdata, params)

    _view_pca(mdata, params)

    _view_umap(mdata, params)

    _view_orgs_by_cluster(mdata, params)

    # If a metadata column was selected
    if params.compare_by:

        _view_grouping_across_clusters(mdata, params)

        _view_orgs_assoc_with_clusters(mdata, params)

        _view_top_assoc_org(mdata, params)
