# THIS FILE IS AUTOGENERATED

from mudata_explorer.sdk.helpers import collapse_params
from mudata_explorer.helpers.add_view import add_view
from muon import MuData


def group_freq(
    mdata: MuData,
    data_table_sidebar=False,
    data_table_axis_value=0,
    data_table_axis_sidebar=False,
    data_table_transforms_value=[],
    data_table_transforms_sidebar=False,
    data_table_tables_value=[],
    data_table_tables_sidebar=False,
    data_table_filter_cols_sidebar=False,
    data_table_filter_cols_type_value=None,
    data_table_filter_cols_type_sidebar=False,
    data_table_filter_cols_tables_value=[],
    data_table_filter_cols_tables_sidebar=False,
    data_table_filter_cols_cname_value=None,
    data_table_filter_cols_cname_sidebar=False,
    data_table_filter_cols_expr_value=None,
    data_table_filter_cols_expr_sidebar=False,
    data_table_filter_cols_value_enum_value=None,
    data_table_filter_cols_value_enum_sidebar=False,
    data_table_filter_cols_value_str_value=None,
    data_table_filter_cols_value_str_sidebar=False,
    data_table_filter_rows_sidebar=False,
    data_table_filter_rows_type_value=None,
    data_table_filter_rows_type_sidebar=False,
    data_table_filter_rows_tables_value=[],
    data_table_filter_rows_tables_sidebar=False,
    data_table_filter_rows_cname_value=None,
    data_table_filter_rows_cname_sidebar=False,
    data_table_filter_rows_expr_value=None,
    data_table_filter_rows_expr_sidebar=False,
    data_table_filter_rows_value_enum_value=None,
    data_table_filter_rows_value_enum_sidebar=False,
    data_table_filter_rows_value_str_value=None,
    data_table_filter_rows_value_str_sidebar=False,
    options_sort_value='Frequency',
    options_sort_sidebar=True,
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
            'data.table.sidebar': extra_params.get('data_table_sidebar', data_table_sidebar),
            'data.table.axis.value': extra_params.get('data_table_axis_value', data_table_axis_value),
            'data.table.axis.sidebar': extra_params.get('data_table_axis_sidebar', data_table_axis_sidebar),
            'data.table.transforms.value': extra_params.get('data_table_transforms_value', data_table_transforms_value),
            'data.table.transforms.sidebar': extra_params.get('data_table_transforms_sidebar', data_table_transforms_sidebar),
            'data.table.tables.value': extra_params.get('data_table_tables_value', data_table_tables_value),
            'data.table.tables.sidebar': extra_params.get('data_table_tables_sidebar', data_table_tables_sidebar),
            'data.table.filter_cols.sidebar': extra_params.get('data_table_filter_cols_sidebar', data_table_filter_cols_sidebar),
            'data.table.filter_cols.type.value': extra_params.get('data_table_filter_cols_type_value', data_table_filter_cols_type_value),
            'data.table.filter_cols.type.sidebar': extra_params.get('data_table_filter_cols_type_sidebar', data_table_filter_cols_type_sidebar),
            'data.table.filter_cols.tables.value': extra_params.get('data_table_filter_cols_tables_value', data_table_filter_cols_tables_value),
            'data.table.filter_cols.tables.sidebar': extra_params.get('data_table_filter_cols_tables_sidebar', data_table_filter_cols_tables_sidebar),
            'data.table.filter_cols.cname.value': extra_params.get('data_table_filter_cols_cname_value', data_table_filter_cols_cname_value),
            'data.table.filter_cols.cname.sidebar': extra_params.get('data_table_filter_cols_cname_sidebar', data_table_filter_cols_cname_sidebar),
            'data.table.filter_cols.expr.value': extra_params.get('data_table_filter_cols_expr_value', data_table_filter_cols_expr_value),
            'data.table.filter_cols.expr.sidebar': extra_params.get('data_table_filter_cols_expr_sidebar', data_table_filter_cols_expr_sidebar),
            'data.table.filter_cols.value_enum.value': extra_params.get('data_table_filter_cols_value_enum_value', data_table_filter_cols_value_enum_value),
            'data.table.filter_cols.value_enum.sidebar': extra_params.get('data_table_filter_cols_value_enum_sidebar', data_table_filter_cols_value_enum_sidebar),
            'data.table.filter_cols.value_str.value': extra_params.get('data_table_filter_cols_value_str_value', data_table_filter_cols_value_str_value),
            'data.table.filter_cols.value_str.sidebar': extra_params.get('data_table_filter_cols_value_str_sidebar', data_table_filter_cols_value_str_sidebar),
            'data.table.filter_rows.sidebar': extra_params.get('data_table_filter_rows_sidebar', data_table_filter_rows_sidebar),
            'data.table.filter_rows.type.value': extra_params.get('data_table_filter_rows_type_value', data_table_filter_rows_type_value),
            'data.table.filter_rows.type.sidebar': extra_params.get('data_table_filter_rows_type_sidebar', data_table_filter_rows_type_sidebar),
            'data.table.filter_rows.tables.value': extra_params.get('data_table_filter_rows_tables_value', data_table_filter_rows_tables_value),
            'data.table.filter_rows.tables.sidebar': extra_params.get('data_table_filter_rows_tables_sidebar', data_table_filter_rows_tables_sidebar),
            'data.table.filter_rows.cname.value': extra_params.get('data_table_filter_rows_cname_value', data_table_filter_rows_cname_value),
            'data.table.filter_rows.cname.sidebar': extra_params.get('data_table_filter_rows_cname_sidebar', data_table_filter_rows_cname_sidebar),
            'data.table.filter_rows.expr.value': extra_params.get('data_table_filter_rows_expr_value', data_table_filter_rows_expr_value),
            'data.table.filter_rows.expr.sidebar': extra_params.get('data_table_filter_rows_expr_sidebar', data_table_filter_rows_expr_sidebar),
            'data.table.filter_rows.value_enum.value': extra_params.get('data_table_filter_rows_value_enum_value', data_table_filter_rows_value_enum_value),
            'data.table.filter_rows.value_enum.sidebar': extra_params.get('data_table_filter_rows_value_enum_sidebar', data_table_filter_rows_value_enum_sidebar),
            'data.table.filter_rows.value_str.value': extra_params.get('data_table_filter_rows_value_str_value', data_table_filter_rows_value_str_value),
            'data.table.filter_rows.value_str.sidebar': extra_params.get('data_table_filter_rows_value_str_sidebar', data_table_filter_rows_value_str_sidebar),
            'options.sort.value': extra_params.get('options_sort_value', options_sort_value),
            'options.sort.sidebar': extra_params.get('options_sort_sidebar', options_sort_sidebar)
        }
    )
