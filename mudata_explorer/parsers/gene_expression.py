from anndata import AnnData
from typing import List, Optional
from mudata import MuData
import streamlit as st
from mudata_explorer.helpers.cirro_readers import util
from mudata_explorer.parsers import common
from mudata_explorer.sdk import view
from dataclasses import dataclass


def parse_adata(adata: AnnData) -> Optional[MuData]:

    mdata = common.parse_adata(
        adata,
        groupby_var=False,
        sum_to_one=False,
        kw="expression"
    )

    params = _get_params(mdata)
    if params is None:
        return

    if st.button("Load Dataset", key="load-microbiome-dataset"):
        _run_processes(mdata, params)
        with st.spinner("Configuring Displays"):
            _add_views(mdata, params)
        return mdata
    

@dataclass
class GeneExpressionParams(common.CompositionalParams):
    """Configuration for the analysis."""
    # Name of the dataset, to be used in the title
    dataset_name: str = "Gene Expression Analysis"
    # Cluster samples with k-means
    k: int = 3
    # Specific genes to plot
    plot_genes: str = ""
    # The comparison may only be categorical
    force_categorical: bool = True
    # Group in the comparison column to be used as a reference
    ref_group: Optional[str] = None
    # Group in the comparison column to be compared to the reference
    comp_group: Optional[str] = None
    # Alpha value for DESeq2
    deseq2_alpha: float = 0.05


def _get_params(mdata: MuData) -> Optional[GeneExpressionParams]:
    return GeneExpressionParams.get_params(mdata)


def _run_processes(mdata: MuData, params: GeneExpressionParams):

    # Get summary metrics for each feature
    with st.spinner("Calculating Gene Summary Statistics"):
        util.summarize_features(mdata, mod="expression")

    # Run PCA
    with st.spinner("Running PCA"):
        util.run_pca(
            mdata,
            mod="expression",
            dest_key="pca"
        )

    # Run kmeans clustering on the samples
    with st.spinner(f"Running K-Means Clustering (k={params.k})"):
        util.run_kmeans(
            mdata,
            mod="expression",
            k=params.k,
            dest_key="kmeans"
        )

    # Run Kruskal to find the features which vary across clusters
    with st.spinner("Finding Genes Driving Clusters"):
        util.run_kruskal(
            mdata,
            table="expression.data",
            dest_key="kruskal_kmeans",
            grouping_cname="kmeans",
            grouping_table="Observation Metadata"
        )

    # If the user selected a metadata column
    if params.compare_by:

        # If a comparison and reference group were selected
        if params.ref_group and params.comp_group:

            # Run DESeq2 to find the features which vary across groups
            with st.spinner(f"Comparing {params.ref_group} and {params.comp_group}"):
                util.run_deseq2(
                    mdata,
                    table="expression.data",
                    dest_key=_deseq2_table_name(params),
                    grouping_cname=params.compare_by,
                    grouping_table="Observation Metadata",
                    ref_level=params.ref_group,
                    comp_level=params.comp_group,
                    alpha=params.deseq2_alpha
                )

        # If no specific comparisons were selected
        else:
            
            # Run Kruskal to find the features which vary across groups
            with st.spinner(f"Comparing all groups of {params.label}"):
                util.run_kruskal(
                    mdata,
                    table="expression.data",
                    dest_key=f"kruskal_{params.compare_by}",
                    grouping_cname=params.compare_by,
                    grouping_table="Observation Metadata"
                )


def _view_pca(mdata: MuData, params: GeneExpressionParams):
    pc1_cname = mdata.mod["expression"].obsm["pca"].columns[0]
    pc2_cname = mdata.mod["expression"].obsm["pca"].columns[1]

    # Get the organisms with the largest loadings on PC1 and PC2
    top_orgs = list(
        mdata.mod["expression"]
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
        data_columns_x_table_value="expression.obsm.pca",
        data_columns_y_table_value="expression.obsm.pca",
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
        vectors_columns_x_table_value="expression.varm.pca",
        vectors_columns_y_table_value="expression.varm.pca",
        vectors_columns_x_cname_value=pc1_cname,
        vectors_columns_x_label_value=pc1_cname,
        vectors_columns_y_cname_value=pc2_cname,
        vectors_columns_y_label_value=pc1_cname,
        vectors_columns_label_enabled_value=False,
        vectors_filter_rows_type_value="index",
        vectors_filter_rows_tables_value=["expression.varm.pca"],
        vectors_filter_rows_expr_value="in",
        vectors_filter_rows_value_enum_value=top_orgs,
        vectors_filter_rows_value_enum_sidebar=True,
        formatting_title_value="Beta Diversity PCA",
        formatting_legend_value=f"""Each point represents a single sample.
The samples are arranged in a two-dimensional space using the PCA algorithm
such that samples with similar patterns of features are closer together.
The color of each point represents the annotated value of '{params.label}'.""",
        formatting_legend_sidebar=True
    )


def _view_kruskal_results(mdata: MuData, params: GeneExpressionParams):
    """
    Summarize the abundances of the genes which are most
    different across the user-specified groups, showing
    their mean abundance in each group.
    """

    view.plotly_category_summarize_values(
        mdata,
        table_category_columns_category_table_value="Observation Metadata",
        table_category_columns_category_cname_value="kmeans",
        table_category_columns_category_label_value="K-Means Cluster",
        table_data_tables_value=["expression.data"],
        table_data_filter_cols=common.top_features_filter_cols(
            mdata,
            params.n_top_features,
            varm_key=f"kruskal_{params.compare_by}",
            cname="rank",
            ascending=True,
            mod="expression"
        ),
        table_data_transforms_value=["zscores_cols"],
        formatting_sort_by_value="Values",
        formatting_color_value="None",
        formatting_title_value=f"Most Variable Genes Across {params.label} Groups",
        formatting_legend_value=f"""The table above shows the top genes which
are the most different between the groups of '{params.label}'.
The expression of each gene is shown for each of those groups.
The size of each point represents the mean expression of the gene among the
samples from that particular group.
To account for different scales of expression, the values have been transformed
to z-scores.""",
        formatting_legend_sidebar=True
    )


def _view_genes_by_cluster(mdata: MuData, params: GeneExpressionParams):

    # Show the most variable organisms by cluster
    view.plotly_category_summarize_values(
        mdata,
        table_category_columns_category_table_value="Observation Metadata",
        table_category_columns_category_cname_value="kmeans",
        table_category_columns_category_label_value="K-Means Cluster",
        table_data_tables_value=["expression.data"],
        table_data_filter_cols=common.top_features_filter_cols(
            mdata,
            params.n_top_features,
            varm_key="summary_stats",
            cname="std",
            ascending=False,
            mod="expression"
        ),
        table_data_transforms_value=["zscores_cols"],
        formatting_sort_by_value="Values",
        formatting_color_value="None",
        formatting_title_value="Most Variable Genes Across Clusters",
        formatting_legend_value="""The table above shows the top genes which
have the highest standard deviation of expression across all samples.
The expression of each gene is shown for each of the clusters identified by the
K-Means algorithm.
The size of each point represents the mean expression of the gene among the
samples which were assigned to a particular cluster.
To account for different scales of expression, the values have been transformed
to z-scores.""",
        formatting_legend_sidebar=True
    )


def _view_grouping_across_clusters(mdata: MuData, params: GeneExpressionParams):

    # Comparison of the variable of interest across clusters
    title = f"Distribution of {params.label} across Clusters"
    legend = f"""The distribution of the variable '{params.label}'
is shown across each of the clusters identified by the K-Means algorithm.
"""

    x_kwargs = dict(
        table="Observation Metadata",
        cname="kmeans",
        label="K-Means Cluster",
    )
    y_kwargs = dict(
        table="Observation Metadata",
        cname=params.compare_by,
        label=params.label
    )

    # For categorical variables, show a category
    # count which includes the k-means clusters
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


def _deseq2_table_name(params: GeneExpressionParams):
    return f"deseq2_{params.compare_by}_{params.ref_group}_{params.comp_group}"


def _view_deseq2_results(mdata: MuData, params: GeneExpressionParams):

    stats_table = _deseq2_table_name(params)

    # Make a volcano plot
    util.add_scatter(
        mdata,
        title=f"Volcano Plot ({params.ref_group} vs {params.comp_group})",
        legend=f"""A summary of DESeq2 results comparing the '{params.ref_group}'
and '{params.comp_group}' groups. The log2 fold change is shown on the x-axis
and the -log10 p-value is shown on the y-axis. The points are colored based on
the mean expression across all samples.""",
        x="log2FoldChange",
        xlabel="Fold Change (log2)",
        y="neg_log10_pvalue",
        ylabel="-log10(p-value)",
        table=f"expression.varm.{stats_table}",
        axis=1,
        color_table=f"expression.varm.{stats_table}",
        cname="baseMean",
        label="Mean Expression",
        is_categorical=False,
        scale="bluered"
    )

    # Make a MA plot
    util.add_scatter(
        mdata,
        title=f"MA Plot ({params.ref_group} vs {params.comp_group})",
        legend=f"""A summary of DESeq2 results comparing the '{params.ref_group}'
and '{params.comp_group}' groups. The mean abundance is shown on the x-axis and
the log2 fold change is shown on the y-axis. The points are colored based on
the -log10(p-value) for the DESeq2 results.""",
        x="baseMean",
        xlabel="Mean Expression",
        y="log2FoldChange",
        ylabel="Fold Change (log2)",
        table=f"expression.varm.{stats_table}",
        axis=1,
        color_table=f"expression.varm.{stats_table}",
        cname="neg_log10_pvalue",
        label="-log10(p-value)",
        is_categorical=False,
        scale="bluered",
        log_x=True
    )

    # Show the table of results
    view.table(
        mdata,
        options_sort_axis_value=1,
        options_sort_columns_sort_by_table_value=f"expression.varm.{stats_table}",
        options_sort_columns_sort_by_cname_value="pvalue",
        options_sort_columns_sort_by_label_value="pvalue",
        data_table_axis_value=1,
        data_table_tables_value=[f"expression.varm.{stats_table}"],
        display_options_title_value=f"DESeq2 Results ({params.ref_group} vs {params.comp_group})"
    )

    # Plot the top result
    _view_specific_gene(
        mdata,
        params,
        (
            mdata.mod["expression"]
            .varm[stats_table]
            .sort_values(by="pvalue")
            .index[0]
        )
    )


def _view_specific_genes(mdata: MuData, params: GeneExpressionParams):
    # For each of the genes
    for gene in params.plot_genes.split(","):
        gene = gene.strip()
        if gene:

            _view_specific_gene(mdata, params, gene)


def _view_specific_gene(mdata: MuData, params: GeneExpressionParams, gene: str):

    # If a comparison column was selected
    if params.compare_by:

        # Show a boxplot with expression across the groups
        view.plotly_box(
            mdata,
            formatting_title_value=f"Expression of {gene} Across {params.label} Groups",
            formatting_legend_value=f"""The boxplot above shows the expression of the gene '{gene}'
across the different groups of '{params.label}'.""",
            formatting_legend_sidebar=True,
            data_axis_value=0,
            data_columns_x_table_value="Observation Metadata",
            data_columns_x_cname_value=params.compare_by,
            data_columns_x_label_value=params.label,
            data_columns_y_table_value="expression.data",
            data_columns_y_cname_value=gene,
            data_columns_y_cname_sidebar=True,
            data_columns_y_label_value=gene,
            data_columns_color_enabled_value=False
        )

    # If no comparison column was selected, show the expression across the k-means clusters
    else:

        # Show a boxplot with expression across the groups
        view.plotly_box(
            mdata,
            formatting_title_value=f"Expression of {gene} Across K-Means Clusters",
            formatting_legend_value=f"""The boxplot above shows the expression of the gene '{gene}'
across the different groups of samples identified by k-means clustering.""",
            formatting_legend_sidebar=True,
            data_axis_value=0,
            data_columns_x_table_value="Observation Metadata",
            data_columns_x_cname_value="kmeans",
            data_columns_x_label_value="K-Means Cluster",
            data_columns_y_table_value="expression.data",
            data_columns_y_cname_value=gene,
            data_columns_y_cname_sidebar=True,
            data_columns_y_label_value=gene,
            data_columns_color_enabled_value=False
        )


def _add_views(mdata: MuData, params: GeneExpressionParams):

    if params.dataset_name:
        view.markdown(mdata, f"### {params.dataset_name}", text_sidebar=True)

    _view_pca(mdata, params)

    _view_genes_by_cluster(mdata, params)

    # If the user selected a metadata column
    if params.compare_by:

        # Show the distribution of the k-means clusters
        # across the different groups
        _view_grouping_across_clusters(mdata, params)

        # If a comparison and reference group were selected
        if params.ref_group and params.comp_group:

            _view_deseq2_results(mdata, params)

        else:

            _view_kruskal_results(mdata, params)

    # Plot any specific genes that the user selected
    if params.plot_genes:
        _view_specific_genes(mdata, params)
