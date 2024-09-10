# THIS FILE IS AUTOGENERATED

from mudata_explorer.sdk.helpers import collapse_params
from mudata_explorer.helpers.add_view import add_view
from mudata import MuData


def plotly_category_count(
    mdata: MuData,
    data_sidebar=False,
    data_axis_value=None,
    data_axis_sidebar=False,
    data_transforms_value=[],
    data_transforms_sidebar=False,
    data_columns_x_sidebar=False,
    data_columns_x_table_value=None,
    data_columns_x_table_sidebar=False,
    data_columns_x_cname_value=None,
    data_columns_x_cname_sidebar=False,
    data_columns_x_label_value='x-axis',
    data_columns_x_label_sidebar=False,
    data_columns_x_scale_value=None,
    data_columns_x_scale_sidebar=False,
    data_columns_x_colorscale=False,
    data_columns_x_is_categorical_value=False,
    data_columns_x_is_categorical_sidebar=False,
    data_columns_color_enabled_value=True,
    data_columns_color_enabled_sidebar=False,
    data_columns_color_sidebar=False,
    data_columns_color_table_value=None,
    data_columns_color_table_sidebar=False,
    data_columns_color_cname_value=None,
    data_columns_color_cname_sidebar=False,
    data_columns_color_label_value='Color',
    data_columns_color_label_sidebar=False,
    data_columns_color_scale_value=None,
    data_columns_color_scale_sidebar=False,
    data_columns_color_colorscale=True,
    data_columns_color_is_categorical_value=True,
    data_columns_color_is_categorical_sidebar=False,
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
    ylabel_value='Count',
    ylabel_sidebar=False,
    barmode_value=None,
    barmode_sidebar=True,
    scale_options_log_y_value=None,
    scale_options_log_y_sidebar=True,
    formatting_title_value='',
    formatting_title_sidebar=True,
    formatting_legend_value=None,
    formatting_legend_sidebar=False,
    annotation_options_show_values_value=None,
    annotation_options_show_values_sidebar=False,
    annotation_options_chisquare_value=None,
    annotation_options_chisquare_sidebar=True,
    **extra_params
):
    """
    
    Show a bar graph which depicts the values in a single column,
    showing the number of times that each unique value is present.
    This is a frequency histogram showing categorical values.

    Optionally, a second column can be used to color the bars.
    When the color column is used, the position may be set to either
    stack vertically, or to group the bars side by side.
    
    """

    assert isinstance(mdata, MuData), "mdata must be a MuData object"
    extra_params = collapse_params(extra_params)

    add_view(
        'plotly-category-count',
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
            'data.columns.color.enabled.value': extra_params.get('data_columns_color_enabled_value', data_columns_color_enabled_value),
            'data.columns.color.enabled.sidebar': extra_params.get('data_columns_color_enabled_sidebar', data_columns_color_enabled_sidebar),
            'data.columns.color.sidebar': extra_params.get('data_columns_color_sidebar', data_columns_color_sidebar),
            'data.columns.color.table.value': extra_params.get('data_columns_color_table_value', data_columns_color_table_value),
            'data.columns.color.table.sidebar': extra_params.get('data_columns_color_table_sidebar', data_columns_color_table_sidebar),
            'data.columns.color.cname.value': extra_params.get('data_columns_color_cname_value', data_columns_color_cname_value),
            'data.columns.color.cname.sidebar': extra_params.get('data_columns_color_cname_sidebar', data_columns_color_cname_sidebar),
            'data.columns.color.label.value': extra_params.get('data_columns_color_label_value', data_columns_color_label_value),
            'data.columns.color.label.sidebar': extra_params.get('data_columns_color_label_sidebar', data_columns_color_label_sidebar),
            'data.columns.color.scale.value': extra_params.get('data_columns_color_scale_value', data_columns_color_scale_value),
            'data.columns.color.scale.sidebar': extra_params.get('data_columns_color_scale_sidebar', data_columns_color_scale_sidebar),
            'data.columns.color.colorscale': extra_params.get('data_columns_color_colorscale', data_columns_color_colorscale),
            'data.columns.color.is_categorical.value': extra_params.get('data_columns_color_is_categorical_value', data_columns_color_is_categorical_value),
            'data.columns.color.is_categorical.sidebar': extra_params.get('data_columns_color_is_categorical_sidebar', data_columns_color_is_categorical_sidebar),
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
            'ylabel.value': extra_params.get('ylabel_value', ylabel_value),
            'ylabel.sidebar': extra_params.get('ylabel_sidebar', ylabel_sidebar),
            'barmode.value': extra_params.get('barmode_value', barmode_value),
            'barmode.sidebar': extra_params.get('barmode_sidebar', barmode_sidebar),
            'scale_options.log_y.value': extra_params.get('scale_options_log_y_value', scale_options_log_y_value),
            'scale_options.log_y.sidebar': extra_params.get('scale_options_log_y_sidebar', scale_options_log_y_sidebar),
            'formatting.title.value': extra_params.get('formatting_title_value', formatting_title_value),
            'formatting.title.sidebar': extra_params.get('formatting_title_sidebar', formatting_title_sidebar),
            'formatting.legend.value': extra_params.get('formatting_legend_value', formatting_legend_value),
            'formatting.legend.sidebar': extra_params.get('formatting_legend_sidebar', formatting_legend_sidebar),
            'annotation_options.show_values.value': extra_params.get('annotation_options_show_values_value', annotation_options_show_values_value),
            'annotation_options.show_values.sidebar': extra_params.get('annotation_options_show_values_sidebar', annotation_options_show_values_sidebar),
            'annotation_options.chisquare.value': extra_params.get('annotation_options_chisquare_value', annotation_options_chisquare_value),
            'annotation_options.chisquare.sidebar': extra_params.get('annotation_options_chisquare_sidebar', annotation_options_chisquare_sidebar)
        }
    )
