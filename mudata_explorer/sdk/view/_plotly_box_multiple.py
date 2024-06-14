# THIS FILE IS AUTOGENERATED

from mudata_explorer import app
from mudata_explorer.sdk.helpers import collapse_params
from muon import MuData


def plotly_box_multiple(
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
    table_category_axis=0,
    table_category_category_table=None,
    table_category_category_cname=None,
    table_category_category_label='Category',
    table_category_rows_query_query_type='',
    table_category_rows_query_query_table='',
    table_category_rows_query_query_cname='',
    table_category_rows_query_query_expr='',
    table_category_rows_query_query_value='',
    table_category_cols_query_query_type='',
    table_category_cols_query_query_table='',
    table_category_cols_query_query_cname='',
    table_category_cols_query_query_expr='',
    table_category_cols_query_query_value='',
    table_category_transforms=[],
    scale_options_log_y=None,
    display_options_ncols=1,
    **extra_params
):
    """
    
    Display multiple columns of data as a box graph using Plotly, summarizing
    the data in terms of the median, quartiles, and outliers.

    A collection of columns are used to define the values on the y-axis, and a
    second column is used for the categorical groups which are
    displayed on the x-axis.
    
    """

    assert isinstance(mdata, MuData), "mdata must be a MuData object"
    extra_params = collapse_params(extra_params)

    app.add_view(
        'plotly-box-multiple',
        mdata,
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
            'table.category.axis': extra_params.get('table_category_axis', table_category_axis),
            'table.category.category.table': extra_params.get('table_category_category_table', table_category_category_table),
            'table.category.category.cname': extra_params.get('table_category_category_cname', table_category_category_cname),
            'table.category.category.label': extra_params.get('table_category_category_label', table_category_category_label),
            'table.category.rows_query.query.type': extra_params.get('table_category_rows_query_query_type', table_category_rows_query_query_type),
            'table.category.rows_query.query.table': extra_params.get('table_category_rows_query_query_table', table_category_rows_query_query_table),
            'table.category.rows_query.query.cname': extra_params.get('table_category_rows_query_query_cname', table_category_rows_query_query_cname),
            'table.category.rows_query.query.expr': extra_params.get('table_category_rows_query_query_expr', table_category_rows_query_query_expr),
            'table.category.rows_query.query.value': extra_params.get('table_category_rows_query_query_value', table_category_rows_query_query_value),
            'table.category.cols_query.query.type': extra_params.get('table_category_cols_query_query_type', table_category_cols_query_query_type),
            'table.category.cols_query.query.table': extra_params.get('table_category_cols_query_query_table', table_category_cols_query_query_table),
            'table.category.cols_query.query.cname': extra_params.get('table_category_cols_query_query_cname', table_category_cols_query_query_cname),
            'table.category.cols_query.query.expr': extra_params.get('table_category_cols_query_query_expr', table_category_cols_query_query_expr),
            'table.category.cols_query.query.value': extra_params.get('table_category_cols_query_query_value', table_category_cols_query_query_value),
            'table.category.transforms': extra_params.get('table_category_transforms', table_category_transforms),
            'scale_options.log_y': extra_params.get('scale_options_log_y', scale_options_log_y),
            'display_options.ncols': extra_params.get('display_options_ncols', display_options_ncols)
        }
    )
