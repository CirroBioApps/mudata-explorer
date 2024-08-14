# THIS FILE IS AUTOGENERATED

from mudata_explorer.sdk.helpers import collapse_params
from mudata_explorer.helpers.assets import make_process
from muon import MuData


def kmeans(
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
    clustering_k=5,
    clustering_min_k=2,
    clustering_max_k=10,
    outputs_dest_key='kmeans',
    **extra_params
):
    """
    
    K-Means clustering is a method of vector quantization that aims to partition
    n observations into k clusters in which each observation belongs to the
    cluster with the nearest mean.

    The user selects the number of clusters (k) to group the samples into.
    If the number of clusters is set to 0, the algorithm will evaluate the
    clustering performance over a range of values of k and use the silhouette
    score to determine the optimal number of clusters.

    Silhouette scores are computed over a range of values of k to evaluate the
    clustering performance. The silhouette score is a measure of how similar
    an object is to its own cluster compared to other clusters.

    - [Wikipedia: K-Means Clustering](https://en.wikipedia.org/wiki/K-means_clustering)
    - [Visual Explanation of K-Means](https://www.naftaliharris.com/blog/visualizing-k-means-clustering/)
    
    """

    assert isinstance(mdata, MuData), "mdata must be a MuData object"
    extra_params = collapse_params(extra_params)

    # Instantiate the process using all of the parameters
    process = make_process(
        'kmeans',
        params={
            'table.data.axis': extra_params.get('table_data_axis', table_data_axis),
            'table.data.tables': extra_params.get('table_data_tables', table_data_tables),
            'table.data.rows_query.query.type': extra_params.get('table_data_rows_query_query_type', table_data_rows_query_query_type),
            'table.data.rows_query.query.table': extra_params.get('table_data_rows_query_query_table', table_data_rows_query_query_table),
            'table.data.rows_query.query.cname': extra_params.get('table_data_rows_query_query_cname', table_data_rows_query_query_cname),
            'table.data.rows_query.query.expr': extra_params.get('table_data_rows_query_query_expr', table_data_rows_query_query_expr),
            'table.data.rows_query.query.value': extra_params.get('table_data_rows_query_query_value', table_data_rows_query_query_value),
            'table.data.cols_query.query.type': extra_params.get('table_data_cols_query_query_type', table_data_cols_query_query_type),
            'table.data.cols_query.query.table': extra_params.get('table_data_cols_query_query_table', table_data_cols_query_query_table),
            'table.data.cols_query.query.cname': extra_params.get('table_data_cols_query_query_cname', table_data_cols_query_query_cname),
            'table.data.cols_query.query.expr': extra_params.get('table_data_cols_query_query_expr', table_data_cols_query_query_expr),
            'table.data.cols_query.query.value': extra_params.get('table_data_cols_query_query_value', table_data_cols_query_query_value),
            'table.data.transforms': extra_params.get('table_data_transforms', table_data_transforms),
            'clustering.k': extra_params.get('clustering_k', clustering_k),
            'clustering.min_k': extra_params.get('clustering_min_k', clustering_min_k),
            'clustering.max_k': extra_params.get('clustering_max_k', clustering_max_k),
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

