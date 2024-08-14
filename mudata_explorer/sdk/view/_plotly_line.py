# THIS FILE IS AUTOGENERATED

from mudata_explorer.sdk.helpers import collapse_params
from mudata_explorer.helpers.add_view import add_view
from muon import MuData


def plotly_line(
    mdata: MuData,
    data_axis=0,
    data_x_table=None,
    data_x_cname=None,
    data_x_label='x-axis',
    data_y_table=None,
    data_y_cname=None,
    data_y_label='y-axis',
    data_sort_by_table=None,
    data_sort_by_cname=None,
    data_sort_by_label='Sort By',
    data_color_table=None,
    data_color_cname=None,
    data_color_label='Color',
    data_color_enabled=True,
    data_color_is_categorical=False,
    data_color_scale='Viridis',
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
    scale_options_log_x=None,
    scale_options_log_y=None,
    scale_options_title='',
    **extra_params
):
    """
    Display a series of data as a line graph using Plotly.
    """

    assert isinstance(mdata, MuData), "mdata must be a MuData object"
    extra_params = collapse_params(extra_params)

    add_view(
        'plotly-line',
        mdata,
        params={
            'data.axis': extra_params.get('data_axis', data_axis),
            'data.x.table': extra_params.get('data_x_table', data_x_table),
            'data.x.cname': extra_params.get('data_x_cname', data_x_cname),
            'data.x.label': extra_params.get('data_x_label', data_x_label),
            'data.y.table': extra_params.get('data_y_table', data_y_table),
            'data.y.cname': extra_params.get('data_y_cname', data_y_cname),
            'data.y.label': extra_params.get('data_y_label', data_y_label),
            'data.sort_by.table': extra_params.get('data_sort_by_table', data_sort_by_table),
            'data.sort_by.cname': extra_params.get('data_sort_by_cname', data_sort_by_cname),
            'data.sort_by.label': extra_params.get('data_sort_by_label', data_sort_by_label),
            'data.color.table': extra_params.get('data_color_table', data_color_table),
            'data.color.cname': extra_params.get('data_color_cname', data_color_cname),
            'data.color.label': extra_params.get('data_color_label', data_color_label),
            'data.color.enabled': extra_params.get('data_color_enabled', data_color_enabled),
            'data.color.is_categorical': extra_params.get('data_color_is_categorical', data_color_is_categorical),
            'data.color.scale': extra_params.get('data_color_scale', data_color_scale),
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
            'scale_options.log_x': extra_params.get('scale_options_log_x', scale_options_log_x),
            'scale_options.log_y': extra_params.get('scale_options_log_y', scale_options_log_y),
            'scale_options.title': extra_params.get('scale_options_title', scale_options_title)
        }
    )
