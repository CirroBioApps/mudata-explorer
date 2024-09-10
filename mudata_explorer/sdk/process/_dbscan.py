# THIS FILE IS AUTOGENERATED

from mudata_explorer.sdk.helpers import collapse_params
from mudata_explorer.helpers.assets import make_process
from mudata import MuData


def dbscan(
    mdata: MuData,
    table_data_sidebar=False,
    table_data_axis_value=None,
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
    clustering_eps_value=0.5,
    clustering_eps_sidebar=False,
    clustering_min_samples_value=5,
    clustering_min_samples_sidebar=False,
    clustering_metric_value=None,
    clustering_metric_sidebar=False,
    outputs_dest_key_value='dbscan',
    outputs_dest_key_sidebar=False,
    **extra_params
):
    """
    
    DBSCAN (Density-Based Spatial Clustering of Applications with Noise)
    is a clustering algorithm that groups together points that are
    closely packed together (points with many nearby neighbors), marking
    as outliers points that lie alone in low-density regions.

    The algorithm works by defining neighborhoods around each point and
    grouping points that are within a certain distance of each other.
    Points that are within the neighborhood of a core point are considered
    part of the same cluster. Points that are within the neighborhood of
    a non-core point but are not core points themselves are considered
    border points. Points that are not within the neighborhood of any
    core points are considered outliers.

    - [Wikipedia: DBSCAN](https://en.wikipedia.org/wiki/DBSCAN)
    - [Explanation of DBSCAN (video)](https://www.youtube.com/watch?v=C3r7tGRe2eI)
    
    """

    assert isinstance(mdata, MuData), "mdata must be a MuData object"
    extra_params = collapse_params(extra_params)

    # Instantiate the process using all of the parameters
    process = make_process(
        'dbscan',
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
            'clustering.eps.value': extra_params.get('clustering_eps_value', clustering_eps_value),
            'clustering.eps.sidebar': extra_params.get('clustering_eps_sidebar', clustering_eps_sidebar),
            'clustering.min_samples.value': extra_params.get('clustering_min_samples_value', clustering_min_samples_value),
            'clustering.min_samples.sidebar': extra_params.get('clustering_min_samples_sidebar', clustering_min_samples_sidebar),
            'clustering.metric.value': extra_params.get('clustering_metric_value', clustering_metric_value),
            'clustering.metric.sidebar': extra_params.get('clustering_metric_sidebar', clustering_metric_sidebar),
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

