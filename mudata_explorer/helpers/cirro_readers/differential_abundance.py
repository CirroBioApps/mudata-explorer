from pandas import DataFrame
import streamlit as st
from cirro import DataPortalDataset
from mudata import MuData
from typing import Optional
from mudata_explorer.helpers.cirro_readers import util
from mudata_explorer.sdk import io, view


def read(dataset: DataPortalDataset) -> Optional[MuData]:
    """
    Read datasets produced by CirroBio/nf-differential-abundance.

    The displays configured for this dataset will be:

    - A UMAP plot of the samples, colored by the variable of interest
    - A table of the top features sorted by p-value
    - A volcano plot with all features
    - For the top N features:
        - A UMAP plot of the samples, colored by the feature
        - A plot showing the association of the feature with the variable
    - A UMAP plot of the samples, colored by leiden clusters
    - A dotplot showing which features are associated with each cluster
    - A plot showing the association of leiden clusters with the variable
    """

    mdata = _read_mdata(dataset)
    if mdata is None:
        return

    params = _get_params(mdata)
    if params is None:
        return

    if st.button("Load Dataset", key="load-from-cirro"):
        _run_processes(mdata, params)
        _add_views(mdata, params)
        return mdata


def _read_mdata(dataset: DataPortalDataset) -> Optional[MuData]:

    # Get the list of files from the dataset
    files = util.list_files(
        dataset,
        pattern="anndata/.*\.h5ad" # noqa
    )

    # If there is only one file
    if len(files) == 1:
        h5ad_file = files[0]
    # If there are no files
    elif len(files) == 0:
        st.error("No files found")
        return
    # If there is more than one file
    else:
        # Pick the file from the dataset
        h5ad_file = util.select_file(
            dataset,
            files=files
        )
        # If no file is selected
        if h5ad_file is None:
            return

    adata = util.read_h5ad(h5ad_file)

    for table_key, axis_prefix in [
        ("X_pca", "PC"),
        ("X_umap", "UMAP "),
        ("X_tsne", "t-SNE "),
    ]:
        adata.obsm[table_key] = DataFrame(
            adata.obsm[table_key],
            index=adata.obs_names,
            columns=[
                f"{axis_prefix}{i+1}"
                for i in range(adata.obsm[table_key].shape[1])
            ]
        )

    return io.build_mdata(
        dict(
            abund=adata
        ),
        obs=adata.obs
    )


def _get_params(mdata: MuData) -> Optional[dict]:

    # Get the name of the variable which was used for differential abundance
    obs_kw, varm_kw = _find_differential_abundance_variable(mdata)
    params = dict(obs_kw=obs_kw, varm_kw=varm_kw)

    # Ask the user what the label should be for that variable
    with st.container(border=1):
        st.write(
            f"""**Differential Abundance Variable**

The variable '{obs_kw}' was used to compare samples.
The user may provide a different label used for this variable in the display.
""")
        params["label"] = st.text_input(
            "Variable Label",
            obs_kw.replace("_", " "),
            help=f"Label used in displays for '{obs_kw}'"
        )

    # Ask if the variable is categorical or continuous
    with st.container(border=1):
        params["is_categorical"] = util.ask_if_categorical(
            obs_kw,
            mdata.obs[obs_kw]
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


def _find_differential_abundance_variable(mdata: MuData) -> Optional[str]:
    """
    Find the variable which was used for differential abundance.
    There should be a key in .varm called {kw}_mavolcanoand
    kw should be in the list of variables in .obs.
    """

    suffix = "_volcano"
    varm_keys = mdata.mod["abund"].varm.keys()
    for kw in varm_keys:
        if kw.endswith(suffix) and kw[:-len(suffix)] in mdata.obs.columns:
            obs_kw = kw[:-len(suffix)]

            if obs_kw in varm_keys:
                return obs_kw, obs_kw
            elif f"mu.{obs_kw}" in varm_keys:
                return obs_kw, f"mu.{obs_kw}"

    raise Exception("No valid .varm[*_volcano] found")


def _run_processes(mdata: MuData, params: dict):

    # Run leiden clustering on the samples
    util.run_leiden(
        mdata,
        mod="abund",
        resolution=params["leiden_res"],
        metric="braycurtis",
        n_neighbors=15,
        dest_key="leiden"
    )

    # Run Kruskal to find the features which vary across clusters
    util.run_kruskal(
        mdata,
        table="abund.data",
        dest_key="kruskal_leiden",
        grouping_cname="leiden",
        grouping_table="Observation Metadata"
    )


def _add_views(mdata: MuData, params: dict):

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
        table="abund.obsm.X_umap",
        axis=0,
        x="UMAP 1",
        xlabel="UMAP 1",
        y="UMAP 2",
        ylabel="UMAP 2",
        color_table="Observation Metadata",
        cname=params["obs_kw"],
        label=params["label"],
        is_categorical=params["is_categorical"],
        scale="D3" if params["is_categorical"] else "bluered",
    )

    # Add the top features table
    util.add_table(
        mdata,
        title="**Differential Abundance Analysis**",
        legend=f"""
The significance of association for each feature with the
{params['label']} variable is shown in the table above.
""",
        table=f"abund.varm.{params['varm_kw']}",
        axis=1,
        sort_by="p_value"
    )

    # Add a volcano plot
    util.add_scatter(
        mdata,
        title=f"Volcano Plot - {params['label']}",
        legend=f"""
The significance of association for each feature with the
{params['label']} variable is shown in the plot above.
The horizontal axis shows the effect size of the association,
while the vertical axis shows the significance of the association.
The color of each point represents the mean abundance of the feature
across all samples in the dataset.
""",
        table=f"abund.varm.{params['varm_kw']}",
        axis=1,
        x="est_coef",
        xlabel="Estimated Coefficient of Association",
        y="neg_log10_pvalue",
        ylabel="p-value (-log10)",
        color_table=f"abund.varm.{params['varm_kw']}",
        cname="mean_abund",
        label="Mean Abundance",
        is_categorical=False,
        scale="bluered"
    )

    # For each of the top features, plot the association with the variable
    # And also show a UMAP colored by that feature
    for feature_id, feature_info in (
        mdata
        .mod["abund"]
        .varm[params["varm_kw"]]
        .sort_values(by="p_value")
        .head(params["n_top_features"])
        .iterrows()
    ):
        title = " ".join([
            f"Distribution of {feature_id} across {params['label']}",
            f"(p-value: {util.format_float(feature_info['p_value'])})"
        ])

        # If the feature of interest is categorical
        if params["is_categorical"]:
            # Show a boxplot across the categories
            util.add_boxplot(
                mdata,
                title=title,
                legend="",
                x_table="Observation Metadata",
                x_cname=params["obs_kw"],
                x_label=params["label"],
                y_table="abund.data",
                y_cname=feature_id,
                y_label=feature_id
            )
        # If the feature of interest is continuous
        else:
            # Show a scatter plot of the feature against the variable
            view.plotly_scatter(
                mdata,
                formatting_title=title,
                data=dict(
                    axis=0,
                    x=dict(
                        table=["Observation Metadata"],
                        cname=params["obs_kw"],
                        label=params["label"]
                    ),
                    y=dict(
                        table=["abund.data"],
                        cname=feature_id,
                        label=feature_id
                    ),
                    size=dict(enabled=False),
                    color=dict(enabled=False)
                )
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
        table="abund.obsm.X_umap",
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
        cname=params["obs_kw"],
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
