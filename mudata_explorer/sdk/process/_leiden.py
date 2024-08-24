# THIS FILE IS AUTOGENERATED

from mudata_explorer.sdk.helpers import collapse_params
from mudata_explorer.helpers.assets import make_process
from muon import MuData


def leiden(
    mdata: MuData,
    table_sidebar=False,
    table_data_sidebar=False,
    table_data_axis=0,
    table_data_axis_sidebar=False,
    table_data_tables=[],
    table_data_tables_sidebar=False,
    table_data_rows_query_query_type='',
    table_data_rows_query_query_sidebar=False,
    table_data_rows_query_query_table='',
    table_data_rows_query_query_cname='',
    table_data_rows_query_query_expr='',
    table_data_rows_query_query_value='',
    table_data_cols_query_query_type='',
    table_data_cols_query_query_sidebar=False,
    table_data_cols_query_query_table='',
    table_data_cols_query_query_cname='',
    table_data_cols_query_query_expr='',
    table_data_cols_query_query_value='',
    table_data_transforms=[],
    table_data_transforms_sidebar=False,
    clustering_sidebar=False,
    clustering_resolution_sidebar=False,
    clustering_resolution=1.0,
    clustering_n_neighbors_sidebar=False,
    clustering_n_neighbors=15,
    clustering_metric_sidebar=False,
    clustering_metric='euclidean',
    outputs_sidebar=False,
    outputs_dest_key_sidebar=False,
    outputs_dest_key='leiden',
    **extra_params
):
    """
    
    Leiden clustering is a graph-based clustering algorithm that is
    particularly well-suited to single-cell RNA-seq data. It is an
    improved version of the Louvain algorithm that is more robust to
    resolution parameters and can find clusters of varying sizes.

    The primary parameter used to tune the size of clusters generated
    by Leiden is the `resolution` parameter. Higher values of `resolution`
    will lead to more and smaller clusters, while lower values will lead
    to fewer and larger clusters.

    - [Traag, Waltman, and van Eck, 2019](https://www.nature.com/articles/s41598-019-41695-z)
    
    """

    assert isinstance(mdata, MuData), "mdata must be a MuData object"
    extra_params = collapse_params(extra_params)

    # Instantiate the process using all of the parameters
    process = make_process(
        'leiden',
        params={
            'table.sidebar': extra_params.get('table_sidebar', table_sidebar),
            'table.data.sidebar': extra_params.get('table_data_sidebar', table_data_sidebar),
            'table.data.axis': extra_params.get('table_data_axis', table_data_axis),
            'table.data.axis.sidebar': extra_params.get('table_data_axis_sidebar', table_data_axis_sidebar),
            'table.data.tables': extra_params.get('table_data_tables', table_data_tables),
            'table.data.tables.sidebar': extra_params.get('table_data_tables_sidebar', table_data_tables_sidebar),
            'table.data.rows_query.query.type': extra_params.get('table_data_rows_query_query_type', table_data_rows_query_query_type),
            'table.data.rows_query.query.sidebar': extra_params.get('table_data_rows_query_query_sidebar', table_data_rows_query_query_sidebar),
            'table.data.rows_query.query.table': extra_params.get('table_data_rows_query_query_table', table_data_rows_query_query_table),
            'table.data.rows_query.query.cname': extra_params.get('table_data_rows_query_query_cname', table_data_rows_query_query_cname),
            'table.data.rows_query.query.expr': extra_params.get('table_data_rows_query_query_expr', table_data_rows_query_query_expr),
            'table.data.rows_query.query.value': extra_params.get('table_data_rows_query_query_value', table_data_rows_query_query_value),
            'table.data.cols_query.query.type': extra_params.get('table_data_cols_query_query_type', table_data_cols_query_query_type),
            'table.data.cols_query.query.sidebar': extra_params.get('table_data_cols_query_query_sidebar', table_data_cols_query_query_sidebar),
            'table.data.cols_query.query.table': extra_params.get('table_data_cols_query_query_table', table_data_cols_query_query_table),
            'table.data.cols_query.query.cname': extra_params.get('table_data_cols_query_query_cname', table_data_cols_query_query_cname),
            'table.data.cols_query.query.expr': extra_params.get('table_data_cols_query_query_expr', table_data_cols_query_query_expr),
            'table.data.cols_query.query.value': extra_params.get('table_data_cols_query_query_value', table_data_cols_query_query_value),
            'table.data.transforms': extra_params.get('table_data_transforms', table_data_transforms),
            'table.data.transforms.sidebar': extra_params.get('table_data_transforms_sidebar', table_data_transforms_sidebar),
            'clustering.sidebar': extra_params.get('clustering_sidebar', clustering_sidebar),
            'clustering.resolution.sidebar': extra_params.get('clustering_resolution_sidebar', clustering_resolution_sidebar),
            'clustering.resolution': extra_params.get('clustering_resolution', clustering_resolution),
            'clustering.n_neighbors.sidebar': extra_params.get('clustering_n_neighbors_sidebar', clustering_n_neighbors_sidebar),
            'clustering.n_neighbors': extra_params.get('clustering_n_neighbors', clustering_n_neighbors),
            'clustering.metric.sidebar': extra_params.get('clustering_metric_sidebar', clustering_metric_sidebar),
            'clustering.metric': extra_params.get('clustering_metric', clustering_metric),
            'outputs.sidebar': extra_params.get('outputs_sidebar', outputs_sidebar),
            'outputs.dest_key.sidebar': extra_params.get('outputs_dest_key_sidebar', outputs_dest_key_sidebar),
            'outputs.dest_key': extra_params.get('outputs_dest_key', outputs_dest_key)
        },
        mdata=mdata,
        params_editable=False
    )

    assert process.params_editable is False, "params_editable must be False"
    assert isinstance(process.mdata, MuData), type(process.mdata)

    # Get the data from the object
    process.get_data()

    # Run the process
    process.execute()

