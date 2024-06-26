# THIS FILE IS AUTOGENERATED

from mudata_explorer import app
from mudata_explorer.sdk.helpers import collapse_params
from mudata_explorer.helpers import make_process
from muon import MuData


def umap(
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
    umap_params_n_neighbors=15,
    umap_params_min_dist=0.1,
    umap_params_metric='correlation',
    umap_params_n_components=2,
    outputs_dest_key='X_umap',
    **extra_params
):
    """
    
    UMAP is a fairly flexible non-linear dimension reduction algorithm.
    It seeks to learn the manifold structure of your data and find a low
    dimensional embedding that preserves the essential topological structure
    of that manifold.

    Compared to something like PCA, UMAP is able to capture more complex
    relationships between groups of data points using a smaller number of
    dimensions.
    Instead of having to look through PC1, PC2, PC3, etc. to understand the
    structure of your data, UMAP can often capture the same information in
    just UMAP1 and UMAP2.

    - [UMAP Documentation](https://umap-learn.readthedocs.io/en/latest/)
    - [Understanding UMAP](https://pair-code.github.io/understanding-umap/)
    - [UMAP in Single-Cell Analysis](https://alleninstitute.org/resource/what-is-a-umap/)
    
    """

    assert isinstance(mdata, MuData), "mdata must be a MuData object"
    extra_params = collapse_params(extra_params)

    # Instantiate the process using all of the parameters
    process = make_process(
        'umap',
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
            'umap_params.n_neighbors': extra_params.get('umap_params_n_neighbors', umap_params_n_neighbors),
            'umap_params.min_dist': extra_params.get('umap_params_min_dist', umap_params_min_dist),
            'umap_params.metric': extra_params.get('umap_params_metric', umap_params_metric),
            'umap_params.n_components': extra_params.get('umap_params_n_components', umap_params_n_components),
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

