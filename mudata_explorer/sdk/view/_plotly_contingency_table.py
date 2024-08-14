# THIS FILE IS AUTOGENERATED

from mudata_explorer.sdk.helpers import collapse_params
from mudata_explorer.helpers.add_view import add_view
from muon import MuData


def plotly_contingency_table(
    mdata: MuData,
    data_axis=0,
    data_x_table=None,
    data_x_cname=None,
    data_x_label='Category A',
    data_y_table=None,
    data_y_cname=None,
    data_y_label='Category B',
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
    formatting_colorscale='blues',
    formatting_values='Number of Items',
    **extra_params
):
    """
    
Compare the number of samples which are found in each combination of
two different columns of categories. This is a frequency table which
shows the number of times that each unique combination of values
is found in the data.

The primary purpose of this display is to identify when there is
a strong correlation between the values in two different columns.

The display can be used to show either:

- The number of items in each combination, or
- The odds ratio of each combination, calculated as the log2 of the
    ratio of the observed frequency to the expected frequency.

    
    """

    assert isinstance(mdata, MuData), "mdata must be a MuData object"
    extra_params = collapse_params(extra_params)

    add_view(
        'plotly-contingency-table',
        mdata,
        params={
            'data.axis': extra_params.get('data_axis', data_axis),
            'data.x.table': extra_params.get('data_x_table', data_x_table),
            'data.x.cname': extra_params.get('data_x_cname', data_x_cname),
            'data.x.label': extra_params.get('data_x_label', data_x_label),
            'data.y.table': extra_params.get('data_y_table', data_y_table),
            'data.y.cname': extra_params.get('data_y_cname', data_y_cname),
            'data.y.label': extra_params.get('data_y_label', data_y_label),
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
            'formatting.colorscale': extra_params.get('formatting_colorscale', formatting_colorscale),
            'formatting.values': extra_params.get('formatting_values', formatting_values)
        }
    )
