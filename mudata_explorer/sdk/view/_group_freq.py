# THIS FILE IS AUTOGENERATED

from mudata_explorer.sdk.helpers import collapse_params
from mudata_explorer.helpers.add_view import add_view
from muon import MuData


def group_freq(
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
    options_sort='Frequency',
    **extra_params
):
    """
    Show a table of the frequency of every combination
    of unique values across any number of selected columns.

This display is most effective for a small number of categorical values
in which the frequency of each unique combination is of interest.

Example input:

| Column A | Column B | Column C |
|----------|----------|----------|
| A        | X        | Y        |
| A        | Y        | Y        |
| B        | X        | Y        |
| B        | X        | Y        |
| B        | Y        | Y        |
| B        | Y        | Y        |

Example output:

| Column A | Column B | Column C | Frequency |
|----------|----------|----------|-----------|
| A        | X        | Y        | 1         |
| A        | Y        | Y        | 1         |
| B        | X        | Y        | 2         |
| B        | Y        | Y        | 2         |

    
    """

    assert isinstance(mdata, MuData), "mdata must be a MuData object"
    extra_params = collapse_params(extra_params)

    add_view(
        'group-freq',
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
            'options.sort': extra_params.get('options_sort', options_sort)
        }
    )
