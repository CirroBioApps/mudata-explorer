# THIS FILE IS AUTOGENERATED

from mudata_explorer.sdk.helpers import collapse_params
from mudata_explorer.helpers.assets import make_process
from muon import MuData


def hdbscan(
    mdata: MuData,
    table_data_sidebar=False,
    table_data_axis_value=0,
    table_data_axis_sidebar=False,
    table_data_transforms_value=[],
    table_data_transforms_sidebar=False,
    table_data_tables_value=[],
    table_data_tables_sidebar=False,
    table_data_filter_cols_sidebar=False,
    table_data_filter_cols_type_value=None,
    table_data_filter_cols_type_sidebar=False,
    table_data_filter_cols_tables_value=[],
    table_data_filter_cols_tables_sidebar=False,
    table_data_filter_cols_cname_value=None,
    table_data_filter_cols_cname_sidebar=False,
    table_data_filter_cols_expr_value=None,
    table_data_filter_cols_expr_sidebar=False,
    table_data_filter_cols_value_enum_value=None,
    table_data_filter_cols_value_enum_sidebar=False,
    table_data_filter_cols_value_str_value=None,
    table_data_filter_cols_value_str_sidebar=False,
    table_data_filter_rows_sidebar=False,
    table_data_filter_rows_type_value=None,
    table_data_filter_rows_type_sidebar=False,
    table_data_filter_rows_tables_value=[],
    table_data_filter_rows_tables_sidebar=False,
    table_data_filter_rows_cname_value=None,
    table_data_filter_rows_cname_sidebar=False,
    table_data_filter_rows_expr_value=None,
    table_data_filter_rows_expr_sidebar=False,
    table_data_filter_rows_value_enum_value=None,
    table_data_filter_rows_value_enum_sidebar=False,
    table_data_filter_rows_value_str_value=None,
    table_data_filter_rows_value_str_sidebar=False,
    clustering_min_cluster_size_value=5,
    clustering_min_cluster_size_sidebar=False,
    clustering_min_samples_value=5,
    clustering_min_samples_sidebar=False,
    clustering_cluster_selection_epsilon_value=0.0,
    clustering_cluster_selection_epsilon_sidebar=False,
    clustering_metric_value='euclidean',
    clustering_metric_sidebar=False,
    clustering_alpha_value=1.0,
    clustering_alpha_sidebar=False,
    clustering_algorithm_value='auto”',
    clustering_algorithm_sidebar=False,
    clustering_leaf_size_value=40,
    clustering_leaf_size_sidebar=False,
    clustering_cluster_selection_method_value='eom',
    clustering_cluster_selection_method_sidebar=False,
    clustering_allow_single_cluster_value=False,
    clustering_allow_single_cluster_sidebar=False,
    outputs_dest_key_value='dbscan',
    outputs_dest_key_sidebar=False,
    **extra_params
):
    """
    
    Cluster data using hierarchical density-based clustering.

    [HDBSCAN - Hierarchical Density-Based Spatial Clustering of Applications with Noise](https://scikit-learn.org/stable/modules/generated/sklearn.cluster.HDBSCAN.html#sklearn.cluster.HDBSCAN).
    Performs DBSCAN over varying epsilon values and integrates the result to find a
    clustering that gives the best stability over epsilon. This allows HDBSCAN to
    find clusters of varying densities (unlike DBSCAN), and be more robust to
    parameter selection. Read more in the [User Guide](https://scikit-learn.org/stable/modules/clustering.html#hdbscan).

    For an example of how to use HDBSCAN, as well as a comparison to DBSCAN, please see the
    [plotting demo](https://scikit-learn.org/stable/auto_examples/cluster/plot_hdbscan.html#sphx-glr-auto-examples-cluster-plot-hdbscan-py).
    
    """

    assert isinstance(mdata, MuData), "mdata must be a MuData object"
    extra_params = collapse_params(extra_params)

    # Instantiate the process using all of the parameters
    process = make_process(
        'hdbscan',
        params={
            'table.data.sidebar': extra_params.get('table_data_sidebar', table_data_sidebar),
            'table.data.axis.value': extra_params.get('table_data_axis_value', table_data_axis_value),
            'table.data.axis.sidebar': extra_params.get('table_data_axis_sidebar', table_data_axis_sidebar),
            'table.data.transforms.value': extra_params.get('table_data_transforms_value', table_data_transforms_value),
            'table.data.transforms.sidebar': extra_params.get('table_data_transforms_sidebar', table_data_transforms_sidebar),
            'table.data.tables.value': extra_params.get('table_data_tables_value', table_data_tables_value),
            'table.data.tables.sidebar': extra_params.get('table_data_tables_sidebar', table_data_tables_sidebar),
            'table.data.filter_cols.sidebar': extra_params.get('table_data_filter_cols_sidebar', table_data_filter_cols_sidebar),
            'table.data.filter_cols.type.value': extra_params.get('table_data_filter_cols_type_value', table_data_filter_cols_type_value),
            'table.data.filter_cols.type.sidebar': extra_params.get('table_data_filter_cols_type_sidebar', table_data_filter_cols_type_sidebar),
            'table.data.filter_cols.tables.value': extra_params.get('table_data_filter_cols_tables_value', table_data_filter_cols_tables_value),
            'table.data.filter_cols.tables.sidebar': extra_params.get('table_data_filter_cols_tables_sidebar', table_data_filter_cols_tables_sidebar),
            'table.data.filter_cols.cname.value': extra_params.get('table_data_filter_cols_cname_value', table_data_filter_cols_cname_value),
            'table.data.filter_cols.cname.sidebar': extra_params.get('table_data_filter_cols_cname_sidebar', table_data_filter_cols_cname_sidebar),
            'table.data.filter_cols.expr.value': extra_params.get('table_data_filter_cols_expr_value', table_data_filter_cols_expr_value),
            'table.data.filter_cols.expr.sidebar': extra_params.get('table_data_filter_cols_expr_sidebar', table_data_filter_cols_expr_sidebar),
            'table.data.filter_cols.value_enum.value': extra_params.get('table_data_filter_cols_value_enum_value', table_data_filter_cols_value_enum_value),
            'table.data.filter_cols.value_enum.sidebar': extra_params.get('table_data_filter_cols_value_enum_sidebar', table_data_filter_cols_value_enum_sidebar),
            'table.data.filter_cols.value_str.value': extra_params.get('table_data_filter_cols_value_str_value', table_data_filter_cols_value_str_value),
            'table.data.filter_cols.value_str.sidebar': extra_params.get('table_data_filter_cols_value_str_sidebar', table_data_filter_cols_value_str_sidebar),
            'table.data.filter_rows.sidebar': extra_params.get('table_data_filter_rows_sidebar', table_data_filter_rows_sidebar),
            'table.data.filter_rows.type.value': extra_params.get('table_data_filter_rows_type_value', table_data_filter_rows_type_value),
            'table.data.filter_rows.type.sidebar': extra_params.get('table_data_filter_rows_type_sidebar', table_data_filter_rows_type_sidebar),
            'table.data.filter_rows.tables.value': extra_params.get('table_data_filter_rows_tables_value', table_data_filter_rows_tables_value),
            'table.data.filter_rows.tables.sidebar': extra_params.get('table_data_filter_rows_tables_sidebar', table_data_filter_rows_tables_sidebar),
            'table.data.filter_rows.cname.value': extra_params.get('table_data_filter_rows_cname_value', table_data_filter_rows_cname_value),
            'table.data.filter_rows.cname.sidebar': extra_params.get('table_data_filter_rows_cname_sidebar', table_data_filter_rows_cname_sidebar),
            'table.data.filter_rows.expr.value': extra_params.get('table_data_filter_rows_expr_value', table_data_filter_rows_expr_value),
            'table.data.filter_rows.expr.sidebar': extra_params.get('table_data_filter_rows_expr_sidebar', table_data_filter_rows_expr_sidebar),
            'table.data.filter_rows.value_enum.value': extra_params.get('table_data_filter_rows_value_enum_value', table_data_filter_rows_value_enum_value),
            'table.data.filter_rows.value_enum.sidebar': extra_params.get('table_data_filter_rows_value_enum_sidebar', table_data_filter_rows_value_enum_sidebar),
            'table.data.filter_rows.value_str.value': extra_params.get('table_data_filter_rows_value_str_value', table_data_filter_rows_value_str_value),
            'table.data.filter_rows.value_str.sidebar': extra_params.get('table_data_filter_rows_value_str_sidebar', table_data_filter_rows_value_str_sidebar),
            'clustering.min_cluster_size.value': extra_params.get('clustering_min_cluster_size_value', clustering_min_cluster_size_value),
            'clustering.min_cluster_size.sidebar': extra_params.get('clustering_min_cluster_size_sidebar', clustering_min_cluster_size_sidebar),
            'clustering.min_samples.value': extra_params.get('clustering_min_samples_value', clustering_min_samples_value),
            'clustering.min_samples.sidebar': extra_params.get('clustering_min_samples_sidebar', clustering_min_samples_sidebar),
            'clustering.cluster_selection_epsilon.value': extra_params.get('clustering_cluster_selection_epsilon_value', clustering_cluster_selection_epsilon_value),
            'clustering.cluster_selection_epsilon.sidebar': extra_params.get('clustering_cluster_selection_epsilon_sidebar', clustering_cluster_selection_epsilon_sidebar),
            'clustering.metric.value': extra_params.get('clustering_metric_value', clustering_metric_value),
            'clustering.metric.sidebar': extra_params.get('clustering_metric_sidebar', clustering_metric_sidebar),
            'clustering.alpha.value': extra_params.get('clustering_alpha_value', clustering_alpha_value),
            'clustering.alpha.sidebar': extra_params.get('clustering_alpha_sidebar', clustering_alpha_sidebar),
            'clustering.algorithm.value': extra_params.get('clustering_algorithm_value', clustering_algorithm_value),
            'clustering.algorithm.sidebar': extra_params.get('clustering_algorithm_sidebar', clustering_algorithm_sidebar),
            'clustering.leaf_size.value': extra_params.get('clustering_leaf_size_value', clustering_leaf_size_value),
            'clustering.leaf_size.sidebar': extra_params.get('clustering_leaf_size_sidebar', clustering_leaf_size_sidebar),
            'clustering.cluster_selection_method.value': extra_params.get('clustering_cluster_selection_method_value', clustering_cluster_selection_method_value),
            'clustering.cluster_selection_method.sidebar': extra_params.get('clustering_cluster_selection_method_sidebar', clustering_cluster_selection_method_sidebar),
            'clustering.allow_single_cluster.value': extra_params.get('clustering_allow_single_cluster_value', clustering_allow_single_cluster_value),
            'clustering.allow_single_cluster.sidebar': extra_params.get('clustering_allow_single_cluster_sidebar', clustering_allow_single_cluster_sidebar),
            'outputs.dest_key.value': extra_params.get('outputs_dest_key_value', outputs_dest_key_value),
            'outputs.dest_key.sidebar': extra_params.get('outputs_dest_key_sidebar', outputs_dest_key_sidebar)
        },
        mdata=mdata,
        params_editable=False
    )

    assert process.params_editable is False, "params_editable must be False"
    assert isinstance(process.mdata, MuData), type(process.mdata)

    # Populate the params for the process
    process.populate_params()

    # Run the process
    process.execute()

