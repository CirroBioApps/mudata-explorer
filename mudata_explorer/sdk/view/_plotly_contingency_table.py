# THIS FILE IS AUTOGENERATED

from mudata_explorer.sdk.helpers import collapse_params
from mudata_explorer.helpers.add_view import add_view
from muon import MuData


def plotly_contingency_table(
    mdata: MuData,
    data_sidebar=False,
    data_axis_value=0,
    data_axis_sidebar=False,
    data_transforms_value=[],
    data_transforms_sidebar=False,
    data_columns_x_sidebar=False,
    data_columns_x_table_value=None,
    data_columns_x_table_sidebar=False,
    data_columns_x_cname_value=None,
    data_columns_x_cname_sidebar=False,
    data_columns_x_label_value='Category A',
    data_columns_x_label_sidebar=False,
    data_columns_x_scale_value=None,
    data_columns_x_scale_sidebar=False,
    data_columns_x_colorscale=False,
    data_columns_x_is_categorical_value=False,
    data_columns_x_is_categorical_sidebar=False,
    data_columns_y_sidebar=False,
    data_columns_y_table_value=None,
    data_columns_y_table_sidebar=False,
    data_columns_y_cname_value=None,
    data_columns_y_cname_sidebar=False,
    data_columns_y_label_value='Category B',
    data_columns_y_label_sidebar=False,
    data_columns_y_scale_value=None,
    data_columns_y_scale_sidebar=False,
    data_columns_y_colorscale=False,
    data_columns_y_is_categorical_value=False,
    data_columns_y_is_categorical_sidebar=False,
    data_filter_rows_sidebar=False,
    data_filter_rows_type_value=None,
    data_filter_rows_type_sidebar=False,
    data_filter_rows_tables_value=[],
    data_filter_rows_tables_sidebar=False,
    data_filter_rows_cname_value=None,
    data_filter_rows_cname_sidebar=False,
    data_filter_rows_expr_value=None,
    data_filter_rows_expr_sidebar=False,
    data_filter_rows_value_enum_value=None,
    data_filter_rows_value_enum_sidebar=False,
    data_filter_rows_value_str_value=None,
    data_filter_rows_value_str_sidebar=False,
    formatting_colorscale_value='blues',
    formatting_colorscale_sidebar=False,
    formatting_values_value='Number of Items',
    formatting_values_sidebar=True,
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
            'data.sidebar': extra_params.get('data_sidebar', data_sidebar),
            'data.axis.value': extra_params.get('data_axis_value', data_axis_value),
            'data.axis.sidebar': extra_params.get('data_axis_sidebar', data_axis_sidebar),
            'data.transforms.value': extra_params.get('data_transforms_value', data_transforms_value),
            'data.transforms.sidebar': extra_params.get('data_transforms_sidebar', data_transforms_sidebar),
            'data.columns.x.sidebar': extra_params.get('data_columns_x_sidebar', data_columns_x_sidebar),
            'data.columns.x.table.value': extra_params.get('data_columns_x_table_value', data_columns_x_table_value),
            'data.columns.x.table.sidebar': extra_params.get('data_columns_x_table_sidebar', data_columns_x_table_sidebar),
            'data.columns.x.cname.value': extra_params.get('data_columns_x_cname_value', data_columns_x_cname_value),
            'data.columns.x.cname.sidebar': extra_params.get('data_columns_x_cname_sidebar', data_columns_x_cname_sidebar),
            'data.columns.x.label.value': extra_params.get('data_columns_x_label_value', data_columns_x_label_value),
            'data.columns.x.label.sidebar': extra_params.get('data_columns_x_label_sidebar', data_columns_x_label_sidebar),
            'data.columns.x.scale.value': extra_params.get('data_columns_x_scale_value', data_columns_x_scale_value),
            'data.columns.x.scale.sidebar': extra_params.get('data_columns_x_scale_sidebar', data_columns_x_scale_sidebar),
            'data.columns.x.colorscale': extra_params.get('data_columns_x_colorscale', data_columns_x_colorscale),
            'data.columns.x.is_categorical.value': extra_params.get('data_columns_x_is_categorical_value', data_columns_x_is_categorical_value),
            'data.columns.x.is_categorical.sidebar': extra_params.get('data_columns_x_is_categorical_sidebar', data_columns_x_is_categorical_sidebar),
            'data.columns.y.sidebar': extra_params.get('data_columns_y_sidebar', data_columns_y_sidebar),
            'data.columns.y.table.value': extra_params.get('data_columns_y_table_value', data_columns_y_table_value),
            'data.columns.y.table.sidebar': extra_params.get('data_columns_y_table_sidebar', data_columns_y_table_sidebar),
            'data.columns.y.cname.value': extra_params.get('data_columns_y_cname_value', data_columns_y_cname_value),
            'data.columns.y.cname.sidebar': extra_params.get('data_columns_y_cname_sidebar', data_columns_y_cname_sidebar),
            'data.columns.y.label.value': extra_params.get('data_columns_y_label_value', data_columns_y_label_value),
            'data.columns.y.label.sidebar': extra_params.get('data_columns_y_label_sidebar', data_columns_y_label_sidebar),
            'data.columns.y.scale.value': extra_params.get('data_columns_y_scale_value', data_columns_y_scale_value),
            'data.columns.y.scale.sidebar': extra_params.get('data_columns_y_scale_sidebar', data_columns_y_scale_sidebar),
            'data.columns.y.colorscale': extra_params.get('data_columns_y_colorscale', data_columns_y_colorscale),
            'data.columns.y.is_categorical.value': extra_params.get('data_columns_y_is_categorical_value', data_columns_y_is_categorical_value),
            'data.columns.y.is_categorical.sidebar': extra_params.get('data_columns_y_is_categorical_sidebar', data_columns_y_is_categorical_sidebar),
            'data.filter_rows.sidebar': extra_params.get('data_filter_rows_sidebar', data_filter_rows_sidebar),
            'data.filter_rows.type.value': extra_params.get('data_filter_rows_type_value', data_filter_rows_type_value),
            'data.filter_rows.type.sidebar': extra_params.get('data_filter_rows_type_sidebar', data_filter_rows_type_sidebar),
            'data.filter_rows.tables.value': extra_params.get('data_filter_rows_tables_value', data_filter_rows_tables_value),
            'data.filter_rows.tables.sidebar': extra_params.get('data_filter_rows_tables_sidebar', data_filter_rows_tables_sidebar),
            'data.filter_rows.cname.value': extra_params.get('data_filter_rows_cname_value', data_filter_rows_cname_value),
            'data.filter_rows.cname.sidebar': extra_params.get('data_filter_rows_cname_sidebar', data_filter_rows_cname_sidebar),
            'data.filter_rows.expr.value': extra_params.get('data_filter_rows_expr_value', data_filter_rows_expr_value),
            'data.filter_rows.expr.sidebar': extra_params.get('data_filter_rows_expr_sidebar', data_filter_rows_expr_sidebar),
            'data.filter_rows.value_enum.value': extra_params.get('data_filter_rows_value_enum_value', data_filter_rows_value_enum_value),
            'data.filter_rows.value_enum.sidebar': extra_params.get('data_filter_rows_value_enum_sidebar', data_filter_rows_value_enum_sidebar),
            'data.filter_rows.value_str.value': extra_params.get('data_filter_rows_value_str_value', data_filter_rows_value_str_value),
            'data.filter_rows.value_str.sidebar': extra_params.get('data_filter_rows_value_str_sidebar', data_filter_rows_value_str_sidebar),
            'formatting.colorscale.value': extra_params.get('formatting_colorscale_value', formatting_colorscale_value),
            'formatting.colorscale.sidebar': extra_params.get('formatting_colorscale_sidebar', formatting_colorscale_sidebar),
            'formatting.values.value': extra_params.get('formatting_values_value', formatting_values_value),
            'formatting.values.sidebar': extra_params.get('formatting_values_sidebar', formatting_values_sidebar)
        }
    )
