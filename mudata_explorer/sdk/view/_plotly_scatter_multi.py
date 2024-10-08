# THIS FILE IS AUTOGENERATED

from mudata_explorer.sdk.helpers import collapse_params
from mudata_explorer.helpers.add_view import add_view
from mudata import MuData


def plotly_scatter_multi(
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
    data_columns_y_sidebar=False,
    data_columns_y_table_value=None,
    data_columns_y_table_sidebar=False,
    data_columns_y_cname_value=None,
    data_columns_y_cname_sidebar=False,
    data_columns_y_label_value='y-axis',
    data_columns_y_label_sidebar=False,
    data_columns_y_scale_value=None,
    data_columns_y_scale_sidebar=False,
    data_columns_y_colorscale=False,
    data_columns_y_is_categorical_value=False,
    data_columns_y_is_categorical_sidebar=False,
    data_columns_size_enabled_value=True,
    data_columns_size_enabled_sidebar=False,
    data_columns_size_sidebar=False,
    data_columns_size_table_value=None,
    data_columns_size_table_sidebar=False,
    data_columns_size_cname_value=None,
    data_columns_size_cname_sidebar=False,
    data_columns_size_label_value='size',
    data_columns_size_label_sidebar=False,
    data_columns_size_scale_value=None,
    data_columns_size_scale_sidebar=False,
    data_columns_size_colorscale=False,
    data_columns_size_is_categorical_value=False,
    data_columns_size_is_categorical_sidebar=False,
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
    colors_sidebar=False,
    colors_axis_value=None,
    colors_axis_sidebar=False,
    colors_transforms_value=[],
    colors_transforms_sidebar=False,
    colors_tables_value=[],
    colors_tables_sidebar=False,
    colors_filter_cols_sidebar=False,
    colors_filter_cols_type_value=None,
    colors_filter_cols_type_sidebar=False,
    colors_filter_cols_tables_value=[],
    colors_filter_cols_tables_sidebar=False,
    colors_filter_cols_cname_value=None,
    colors_filter_cols_cname_sidebar=False,
    colors_filter_cols_expr_value=None,
    colors_filter_cols_expr_sidebar=False,
    colors_filter_cols_value_enum_value=None,
    colors_filter_cols_value_enum_sidebar=False,
    colors_filter_cols_value_str_value=None,
    colors_filter_cols_value_str_sidebar=False,
    colors_filter_rows_sidebar=False,
    colors_filter_rows_type_value=None,
    colors_filter_rows_type_sidebar=False,
    colors_filter_rows_tables_value=[],
    colors_filter_rows_tables_sidebar=False,
    colors_filter_rows_cname_value=None,
    colors_filter_rows_cname_sidebar=False,
    colors_filter_rows_expr_value=None,
    colors_filter_rows_expr_sidebar=False,
    colors_filter_rows_value_enum_value=None,
    colors_filter_rows_value_enum_sidebar=False,
    colors_filter_rows_value_str_value=None,
    colors_filter_rows_value_str_sidebar=False,
    scale_options_log_x_value=None,
    scale_options_log_x_sidebar=True,
    scale_options_log_y_value=None,
    scale_options_log_y_sidebar=True,
    scale_options_log_color_value=None,
    scale_options_log_color_sidebar=True,
    formatting_colorscale_value='bluered',
    formatting_colorscale_sidebar=False,
    formatting_opacity_value=1.0,
    formatting_opacity_sidebar=True,
    formatting_ncols_value=1,
    formatting_ncols_sidebar=False,
    formatting_color_label_value='Abundance',
    formatting_color_label_sidebar=False,
    formatting_title_value='',
    formatting_title_sidebar=True,
    formatting_legend_value=None,
    formatting_legend_sidebar=False,
    **extra_params
):
    """
    
Display a two dimensional distribution of data using Plotly.
Display the same set of points in multiple panels, each of which
is colored by a different column of data.

In addition to the data used for the x- and y-axes, the size
of the points can be set to represent additional dimensions of the data.

Additional formatting options include setting the opacity of the points
and using a log scale for the x- and y-axes.

    """

    assert isinstance(mdata, MuData), "mdata must be a MuData object"
    extra_params = collapse_params(extra_params)

    add_view(
        'plotly-scatter-multi',
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
            'data.columns.size.enabled.value': extra_params.get('data_columns_size_enabled_value', data_columns_size_enabled_value),
            'data.columns.size.enabled.sidebar': extra_params.get('data_columns_size_enabled_sidebar', data_columns_size_enabled_sidebar),
            'data.columns.size.sidebar': extra_params.get('data_columns_size_sidebar', data_columns_size_sidebar),
            'data.columns.size.table.value': extra_params.get('data_columns_size_table_value', data_columns_size_table_value),
            'data.columns.size.table.sidebar': extra_params.get('data_columns_size_table_sidebar', data_columns_size_table_sidebar),
            'data.columns.size.cname.value': extra_params.get('data_columns_size_cname_value', data_columns_size_cname_value),
            'data.columns.size.cname.sidebar': extra_params.get('data_columns_size_cname_sidebar', data_columns_size_cname_sidebar),
            'data.columns.size.label.value': extra_params.get('data_columns_size_label_value', data_columns_size_label_value),
            'data.columns.size.label.sidebar': extra_params.get('data_columns_size_label_sidebar', data_columns_size_label_sidebar),
            'data.columns.size.scale.value': extra_params.get('data_columns_size_scale_value', data_columns_size_scale_value),
            'data.columns.size.scale.sidebar': extra_params.get('data_columns_size_scale_sidebar', data_columns_size_scale_sidebar),
            'data.columns.size.colorscale': extra_params.get('data_columns_size_colorscale', data_columns_size_colorscale),
            'data.columns.size.is_categorical.value': extra_params.get('data_columns_size_is_categorical_value', data_columns_size_is_categorical_value),
            'data.columns.size.is_categorical.sidebar': extra_params.get('data_columns_size_is_categorical_sidebar', data_columns_size_is_categorical_sidebar),
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
            'colors.sidebar': extra_params.get('colors_sidebar', colors_sidebar),
            'colors.axis.value': extra_params.get('colors_axis_value', colors_axis_value),
            'colors.axis.sidebar': extra_params.get('colors_axis_sidebar', colors_axis_sidebar),
            'colors.transforms.value': extra_params.get('colors_transforms_value', colors_transforms_value),
            'colors.transforms.sidebar': extra_params.get('colors_transforms_sidebar', colors_transforms_sidebar),
            'colors.tables.value': extra_params.get('colors_tables_value', colors_tables_value),
            'colors.tables.sidebar': extra_params.get('colors_tables_sidebar', colors_tables_sidebar),
            'colors.filter_cols.sidebar': extra_params.get('colors_filter_cols_sidebar', colors_filter_cols_sidebar),
            'colors.filter_cols.type.value': extra_params.get('colors_filter_cols_type_value', colors_filter_cols_type_value),
            'colors.filter_cols.type.sidebar': extra_params.get('colors_filter_cols_type_sidebar', colors_filter_cols_type_sidebar),
            'colors.filter_cols.tables.value': extra_params.get('colors_filter_cols_tables_value', colors_filter_cols_tables_value),
            'colors.filter_cols.tables.sidebar': extra_params.get('colors_filter_cols_tables_sidebar', colors_filter_cols_tables_sidebar),
            'colors.filter_cols.cname.value': extra_params.get('colors_filter_cols_cname_value', colors_filter_cols_cname_value),
            'colors.filter_cols.cname.sidebar': extra_params.get('colors_filter_cols_cname_sidebar', colors_filter_cols_cname_sidebar),
            'colors.filter_cols.expr.value': extra_params.get('colors_filter_cols_expr_value', colors_filter_cols_expr_value),
            'colors.filter_cols.expr.sidebar': extra_params.get('colors_filter_cols_expr_sidebar', colors_filter_cols_expr_sidebar),
            'colors.filter_cols.value_enum.value': extra_params.get('colors_filter_cols_value_enum_value', colors_filter_cols_value_enum_value),
            'colors.filter_cols.value_enum.sidebar': extra_params.get('colors_filter_cols_value_enum_sidebar', colors_filter_cols_value_enum_sidebar),
            'colors.filter_cols.value_str.value': extra_params.get('colors_filter_cols_value_str_value', colors_filter_cols_value_str_value),
            'colors.filter_cols.value_str.sidebar': extra_params.get('colors_filter_cols_value_str_sidebar', colors_filter_cols_value_str_sidebar),
            'colors.filter_rows.sidebar': extra_params.get('colors_filter_rows_sidebar', colors_filter_rows_sidebar),
            'colors.filter_rows.type.value': extra_params.get('colors_filter_rows_type_value', colors_filter_rows_type_value),
            'colors.filter_rows.type.sidebar': extra_params.get('colors_filter_rows_type_sidebar', colors_filter_rows_type_sidebar),
            'colors.filter_rows.tables.value': extra_params.get('colors_filter_rows_tables_value', colors_filter_rows_tables_value),
            'colors.filter_rows.tables.sidebar': extra_params.get('colors_filter_rows_tables_sidebar', colors_filter_rows_tables_sidebar),
            'colors.filter_rows.cname.value': extra_params.get('colors_filter_rows_cname_value', colors_filter_rows_cname_value),
            'colors.filter_rows.cname.sidebar': extra_params.get('colors_filter_rows_cname_sidebar', colors_filter_rows_cname_sidebar),
            'colors.filter_rows.expr.value': extra_params.get('colors_filter_rows_expr_value', colors_filter_rows_expr_value),
            'colors.filter_rows.expr.sidebar': extra_params.get('colors_filter_rows_expr_sidebar', colors_filter_rows_expr_sidebar),
            'colors.filter_rows.value_enum.value': extra_params.get('colors_filter_rows_value_enum_value', colors_filter_rows_value_enum_value),
            'colors.filter_rows.value_enum.sidebar': extra_params.get('colors_filter_rows_value_enum_sidebar', colors_filter_rows_value_enum_sidebar),
            'colors.filter_rows.value_str.value': extra_params.get('colors_filter_rows_value_str_value', colors_filter_rows_value_str_value),
            'colors.filter_rows.value_str.sidebar': extra_params.get('colors_filter_rows_value_str_sidebar', colors_filter_rows_value_str_sidebar),
            'scale_options.log_x.value': extra_params.get('scale_options_log_x_value', scale_options_log_x_value),
            'scale_options.log_x.sidebar': extra_params.get('scale_options_log_x_sidebar', scale_options_log_x_sidebar),
            'scale_options.log_y.value': extra_params.get('scale_options_log_y_value', scale_options_log_y_value),
            'scale_options.log_y.sidebar': extra_params.get('scale_options_log_y_sidebar', scale_options_log_y_sidebar),
            'scale_options.log_color.value': extra_params.get('scale_options_log_color_value', scale_options_log_color_value),
            'scale_options.log_color.sidebar': extra_params.get('scale_options_log_color_sidebar', scale_options_log_color_sidebar),
            'formatting.colorscale.value': extra_params.get('formatting_colorscale_value', formatting_colorscale_value),
            'formatting.colorscale.sidebar': extra_params.get('formatting_colorscale_sidebar', formatting_colorscale_sidebar),
            'formatting.opacity.value': extra_params.get('formatting_opacity_value', formatting_opacity_value),
            'formatting.opacity.sidebar': extra_params.get('formatting_opacity_sidebar', formatting_opacity_sidebar),
            'formatting.ncols.value': extra_params.get('formatting_ncols_value', formatting_ncols_value),
            'formatting.ncols.sidebar': extra_params.get('formatting_ncols_sidebar', formatting_ncols_sidebar),
            'formatting.color_label.value': extra_params.get('formatting_color_label_value', formatting_color_label_value),
            'formatting.color_label.sidebar': extra_params.get('formatting_color_label_sidebar', formatting_color_label_sidebar),
            'formatting.title.value': extra_params.get('formatting_title_value', formatting_title_value),
            'formatting.title.sidebar': extra_params.get('formatting_title_sidebar', formatting_title_sidebar),
            'formatting.legend.value': extra_params.get('formatting_legend_value', formatting_legend_value),
            'formatting.legend.sidebar': extra_params.get('formatting_legend_sidebar', formatting_legend_sidebar)
        }
    )
