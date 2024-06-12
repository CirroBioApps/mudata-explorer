# THIS FILE IS AUTOGENERATED

from mudata_explorer import app
from muon import MuData


def dbscan(
    mdata: MuData,
    table_data_axis=0,
    table_data_tables=[],
    table_data_rows_query_query_type='',
    table_data_rows_query_query_table='',
    table_data_rows_query_query_cname='',
    table_data_rows_query_query_expr='',
    table_data_rows_query_query_value='',
    table_data_cols_query_query_type='',
    table_data_cols_query_query_table='',
    table_data_cols_query_query_cname='',
    table_data_cols_query_query_expr='',
    table_data_cols_query_query_value='',
    table_data_transforms=[],
    clustering_use_zscore=None,
    clustering_eps=0.5,
    clustering_min_samples=5,
    clustering_metric=None,
    outputs_dest_key='dbscan'
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

    app.add_process(
        'dbscan',
        mdata,
        params={
            'table.data.axis': table_data_axis,
            'table.data.tables': table_data_tables,
            'table.data.rows_query.query.type': table_data_rows_query_query_type,
            'table.data.rows_query.query.table': table_data_rows_query_query_table,
            'table.data.rows_query.query.cname': table_data_rows_query_query_cname,
            'table.data.rows_query.query.expr': table_data_rows_query_query_expr,
            'table.data.rows_query.query.value': table_data_rows_query_query_value,
            'table.data.cols_query.query.type': table_data_cols_query_query_type,
            'table.data.cols_query.query.table': table_data_cols_query_query_table,
            'table.data.cols_query.query.cname': table_data_cols_query_query_cname,
            'table.data.cols_query.query.expr': table_data_cols_query_query_expr,
            'table.data.cols_query.query.value': table_data_cols_query_query_value,
            'table.data.transforms': table_data_transforms,
            'clustering.use_zscore': clustering_use_zscore,
            'clustering.eps': clustering_eps,
            'clustering.min_samples': clustering_min_samples,
            'clustering.metric': clustering_metric,
            'outputs.dest_key': outputs_dest_key
        }
    )
