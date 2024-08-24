# THIS FILE IS AUTOGENERATED

from mudata_explorer.sdk.helpers import collapse_params
from mudata_explorer.helpers.add_view import add_view
from muon import MuData


def table(
    mdata: MuData,
    data_sidebar=False,
    data_table_sidebar=False,
    data_table_axis=0,
    data_table_axis_sidebar=False,
    data_table_tables=[],
    data_table_tables_sidebar=False,
    data_table_rows_query_query_type='',
    data_table_rows_query_query_sidebar=False,
    data_table_rows_query_query_table='',
    data_table_rows_query_query_cname='',
    data_table_rows_query_query_expr='',
    data_table_rows_query_query_value='',
    data_table_cols_query_query_type='',
    data_table_cols_query_query_sidebar=False,
    data_table_cols_query_query_table='',
    data_table_cols_query_query_cname='',
    data_table_cols_query_query_expr='',
    data_table_cols_query_query_value='',
    data_table_transforms=[],
    data_table_transforms_sidebar=False,
    options_sidebar=False,
    options_sort_sidebar=False,
    options_sort_axis=0,
    options_sort_axis_sidebar=False,
    options_sort_sort_by_table=None,
    options_sort_sort_by_table_sidebar=False,
    options_sort_sort_by_cname=None,
    options_sort_sort_by_cname_sidebar=False,
    options_sort_sort_by_label='sort_by',
    options_sort_sort_by_label_sidebar=False,
    options_sort_rows_query_query_type='',
    options_sort_rows_query_query_sidebar=False,
    options_sort_rows_query_query_table='',
    options_sort_rows_query_query_cname='',
    options_sort_rows_query_query_expr='',
    options_sort_rows_query_query_value='',
    options_sort_cols_query_query_type='',
    options_sort_cols_query_query_sidebar=False,
    options_sort_cols_query_query_table='',
    options_sort_cols_query_query_cname='',
    options_sort_cols_query_query_expr='',
    options_sort_cols_query_query_value='',
    options_sort_transforms=[],
    options_sort_transforms_sidebar=False,
    **extra_params
):
    """
    Show a table of data.
    """

    assert isinstance(mdata, MuData), "mdata must be a MuData object"
    extra_params = collapse_params(extra_params)

    add_view(
        'table',
        mdata,
        params={
            'data.sidebar': extra_params.get('data_sidebar', data_sidebar),
            'data.table.sidebar': extra_params.get('data_table_sidebar', data_table_sidebar),
            'data.table.axis': extra_params.get('data_table_axis', data_table_axis),
            'data.table.axis.sidebar': extra_params.get('data_table_axis_sidebar', data_table_axis_sidebar),
            'data.table.tables': extra_params.get('data_table_tables', data_table_tables),
            'data.table.tables.sidebar': extra_params.get('data_table_tables_sidebar', data_table_tables_sidebar),
            'data.table.rows_query.query.type': extra_params.get('data_table_rows_query_query_type', data_table_rows_query_query_type),
            'data.table.rows_query.query.sidebar': extra_params.get('data_table_rows_query_query_sidebar', data_table_rows_query_query_sidebar),
            'data.table.rows_query.query.table': extra_params.get('data_table_rows_query_query_table', data_table_rows_query_query_table),
            'data.table.rows_query.query.cname': extra_params.get('data_table_rows_query_query_cname', data_table_rows_query_query_cname),
            'data.table.rows_query.query.expr': extra_params.get('data_table_rows_query_query_expr', data_table_rows_query_query_expr),
            'data.table.rows_query.query.value': extra_params.get('data_table_rows_query_query_value', data_table_rows_query_query_value),
            'data.table.cols_query.query.type': extra_params.get('data_table_cols_query_query_type', data_table_cols_query_query_type),
            'data.table.cols_query.query.sidebar': extra_params.get('data_table_cols_query_query_sidebar', data_table_cols_query_query_sidebar),
            'data.table.cols_query.query.table': extra_params.get('data_table_cols_query_query_table', data_table_cols_query_query_table),
            'data.table.cols_query.query.cname': extra_params.get('data_table_cols_query_query_cname', data_table_cols_query_query_cname),
            'data.table.cols_query.query.expr': extra_params.get('data_table_cols_query_query_expr', data_table_cols_query_query_expr),
            'data.table.cols_query.query.value': extra_params.get('data_table_cols_query_query_value', data_table_cols_query_query_value),
            'data.table.transforms': extra_params.get('data_table_transforms', data_table_transforms),
            'data.table.transforms.sidebar': extra_params.get('data_table_transforms_sidebar', data_table_transforms_sidebar),
            'options.sidebar': extra_params.get('options_sidebar', options_sidebar),
            'options.sort.sidebar': extra_params.get('options_sort_sidebar', options_sort_sidebar),
            'options.sort.axis': extra_params.get('options_sort_axis', options_sort_axis),
            'options.sort.axis.sidebar': extra_params.get('options_sort_axis_sidebar', options_sort_axis_sidebar),
            'options.sort.sort_by.table': extra_params.get('options_sort_sort_by_table', options_sort_sort_by_table),
            'options.sort.sort_by.table.sidebar': extra_params.get('options_sort_sort_by_table_sidebar', options_sort_sort_by_table_sidebar),
            'options.sort.sort_by.cname': extra_params.get('options_sort_sort_by_cname', options_sort_sort_by_cname),
            'options.sort.sort_by.cname.sidebar': extra_params.get('options_sort_sort_by_cname_sidebar', options_sort_sort_by_cname_sidebar),
            'options.sort.sort_by.label': extra_params.get('options_sort_sort_by_label', options_sort_sort_by_label),
            'options.sort.sort_by.label.sidebar': extra_params.get('options_sort_sort_by_label_sidebar', options_sort_sort_by_label_sidebar),
            'options.sort.rows_query.query.type': extra_params.get('options_sort_rows_query_query_type', options_sort_rows_query_query_type),
            'options.sort.rows_query.query.sidebar': extra_params.get('options_sort_rows_query_query_sidebar', options_sort_rows_query_query_sidebar),
            'options.sort.rows_query.query.table': extra_params.get('options_sort_rows_query_query_table', options_sort_rows_query_query_table),
            'options.sort.rows_query.query.cname': extra_params.get('options_sort_rows_query_query_cname', options_sort_rows_query_query_cname),
            'options.sort.rows_query.query.expr': extra_params.get('options_sort_rows_query_query_expr', options_sort_rows_query_query_expr),
            'options.sort.rows_query.query.value': extra_params.get('options_sort_rows_query_query_value', options_sort_rows_query_query_value),
            'options.sort.cols_query.query.type': extra_params.get('options_sort_cols_query_query_type', options_sort_cols_query_query_type),
            'options.sort.cols_query.query.sidebar': extra_params.get('options_sort_cols_query_query_sidebar', options_sort_cols_query_query_sidebar),
            'options.sort.cols_query.query.table': extra_params.get('options_sort_cols_query_query_table', options_sort_cols_query_query_table),
            'options.sort.cols_query.query.cname': extra_params.get('options_sort_cols_query_query_cname', options_sort_cols_query_query_cname),
            'options.sort.cols_query.query.expr': extra_params.get('options_sort_cols_query_query_expr', options_sort_cols_query_query_expr),
            'options.sort.cols_query.query.value': extra_params.get('options_sort_cols_query_query_value', options_sort_cols_query_query_value),
            'options.sort.transforms': extra_params.get('options_sort_transforms', options_sort_transforms),
            'options.sort.transforms.sidebar': extra_params.get('options_sort_transforms_sidebar', options_sort_transforms_sidebar)
        }
    )
