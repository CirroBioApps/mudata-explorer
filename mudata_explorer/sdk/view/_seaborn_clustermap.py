# THIS FILE IS AUTOGENERATED

from mudata_explorer import app
from mudata_explorer.sdk.helpers import collapse_params
from muon import MuData


def seaborn_clustermap(
    mdata: MuData,
    data_axis=0,
    data_tables=[],
    data_rows_query_query_type='',
    data_rows_query_query_table='',
    data_rows_query_query_cname='',
    data_rows_query_query_expr='',
    data_rows_query_query_value='',
    data_cols_query_query_type='',
    data_cols_query_query_table='',
    data_cols_query_query_cname='',
    data_cols_query_query_expr='',
    data_cols_query_query_value='',
    data_transforms=[],
    z_score='None',
    **extra_params
):
    """
    Display a heatmap with clustered rows and columnsusing Seaborn.
    """

    assert isinstance(mdata, MuData), "mdata must be a MuData object"
    extra_params = collapse_params(extra_params)

    app.add_view(
        'seaborn-clustermap',
        mdata,
        params={
            'data.axis': extra_params.get('data_axis', data_axis),
            'data.tables': extra_params.get('data_tables', data_tables),
            'data.rows_query.query.type': extra_params.get('data_rows_query_query_type', data_rows_query_query_type),
            'data.rows_query.query.table': extra_params.get('data_rows_query_query_table', data_rows_query_query_table),
            'data.rows_query.query.cname': extra_params.get('data_rows_query_query_cname', data_rows_query_query_cname),
            'data.rows_query.query.expr': extra_params.get('data_rows_query_query_expr', data_rows_query_query_expr),
            'data.rows_query.query.value': extra_params.get('data_rows_query_query_value', data_rows_query_query_value),
            'data.cols_query.query.type': extra_params.get('data_cols_query_query_type', data_cols_query_query_type),
            'data.cols_query.query.table': extra_params.get('data_cols_query_query_table', data_cols_query_query_table),
            'data.cols_query.query.cname': extra_params.get('data_cols_query_query_cname', data_cols_query_query_cname),
            'data.cols_query.query.expr': extra_params.get('data_cols_query_query_expr', data_cols_query_query_expr),
            'data.cols_query.query.value': extra_params.get('data_cols_query_query_value', data_cols_query_query_value),
            'data.transforms': extra_params.get('data_transforms', data_transforms),
            'z_score': extra_params.get('z_score', z_score)
        }
    )
