# THIS FILE IS AUTOGENERATED

from mudata_explorer.sdk.helpers import collapse_params
from mudata_explorer.helpers.add_view import add_view
from muon import MuData


def plotly_box_multiple(
    mdata: MuData,
    table_data_sidebar=False,
    table_data_axis_value=None,
    table_data_axis_sidebar=False,
    table_data_transforms_value=[],
    table_data_transforms_sidebar=False,
    table_data_tables_value=[],
    table_data_tables_sidebar=False,
    table_data_filter_cols_sidebar=False,
    table_data_filter_cols_type_value=None,
    table_data_filter_cols_type_sidebar=False,
    table_data_filter_cols_tables_value=[],
    table_data_filter_cols_tables_sidebar=False,
    table_data_filter_cols_cname_value=None,
    table_data_filter_cols_cname_sidebar=False,
    table_data_filter_cols_expr_value=None,
    table_data_filter_cols_expr_sidebar=False,
    table_data_filter_cols_value_enum_value=None,
    table_data_filter_cols_value_enum_sidebar=False,
    table_data_filter_cols_value_str_value=None,
    table_data_filter_cols_value_str_sidebar=False,
    table_data_filter_rows_sidebar=False,
    table_data_filter_rows_type_value=None,
    table_data_filter_rows_type_sidebar=False,
    table_data_filter_rows_tables_value=[],
    table_data_filter_rows_tables_sidebar=False,
    table_data_filter_rows_cname_value=None,
    table_data_filter_rows_cname_sidebar=False,
    table_data_filter_rows_expr_value=None,
    table_data_filter_rows_expr_sidebar=False,
    table_data_filter_rows_value_enum_value=None,
    table_data_filter_rows_value_enum_sidebar=False,
    table_data_filter_rows_value_str_value=None,
    table_data_filter_rows_value_str_sidebar=False,
    table_category_enabled_value=True,
    table_category_enabled_sidebar=False,
    table_category_sidebar=False,
    table_category_axis_value=None,
    table_category_axis_sidebar=False,
    table_category_transforms_value=[],
    table_category_transforms_sidebar=False,
    table_category_columns_category_sidebar=False,
    table_category_columns_category_table_value=None,
    table_category_columns_category_table_sidebar=False,
    table_category_columns_category_cname_value=None,
    table_category_columns_category_cname_sidebar=False,
    table_category_columns_category_label_value='Category',
    table_category_columns_category_label_sidebar=False,
    table_category_columns_category_scale_value=None,
    table_category_columns_category_scale_sidebar=False,
    table_category_columns_category_colorscale=False,
    table_category_columns_category_is_categorical_value=False,
    table_category_columns_category_is_categorical_sidebar=False,
    table_category_filter_rows_sidebar=False,
    table_category_filter_rows_type_value=None,
    table_category_filter_rows_type_sidebar=False,
    table_category_filter_rows_tables_value=[],
    table_category_filter_rows_tables_sidebar=False,
    table_category_filter_rows_cname_value=None,
    table_category_filter_rows_cname_sidebar=False,
    table_category_filter_rows_expr_value=None,
    table_category_filter_rows_expr_sidebar=False,
    table_category_filter_rows_value_enum_value=None,
    table_category_filter_rows_value_enum_sidebar=False,
    table_category_filter_rows_value_str_value=None,
    table_category_filter_rows_value_str_sidebar=False,
    variable_options_axis_value='Y-Axis',
    variable_options_axis_sidebar=True,
    variable_options_log_values_value=None,
    variable_options_log_values_sidebar=True,
    variable_options_sort_by_value='Mean',
    variable_options_sort_by_sidebar=True,
    category_options_axis_value='Axis',
    category_options_axis_sidebar=False,
    category_options_sort_by_value='Mean',
    category_options_sort_by_sidebar=False,
    display_options_ncols_value=1,
    display_options_ncols_sidebar=False,
    display_options_outliers_enabled_value=True,
    display_options_outliers_enabled_sidebar=True,
    display_options_title_value='',
    display_options_title_sidebar=True,
    display_options_var_label_value='Variable',
    display_options_var_label_sidebar=False,
    display_options_val_label_value='Value',
    display_options_val_label_sidebar=False,
    display_options_height_value=500,
    display_options_height_sidebar=False,
    display_options_legend_value=None,
    display_options_legend_sidebar=False,
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

    add_view(
        'plotly-box-multiple',
        mdata,
        params={
            'table.data.sidebar': extra_params.get('table_data_sidebar', table_data_sidebar),
            'table.data.axis.value': extra_params.get('table_data_axis_value', table_data_axis_value),
            'table.data.axis.sidebar': extra_params.get('table_data_axis_sidebar', table_data_axis_sidebar),
            'table.data.transforms.value': extra_params.get('table_data_transforms_value', table_data_transforms_value),
            'table.data.transforms.sidebar': extra_params.get('table_data_transforms_sidebar', table_data_transforms_sidebar),
            'table.data.tables.value': extra_params.get('table_data_tables_value', table_data_tables_value),
            'table.data.tables.sidebar': extra_params.get('table_data_tables_sidebar', table_data_tables_sidebar),
            'table.data.filter_cols.sidebar': extra_params.get('table_data_filter_cols_sidebar', table_data_filter_cols_sidebar),
            'table.data.filter_cols.type.value': extra_params.get('table_data_filter_cols_type_value', table_data_filter_cols_type_value),
            'table.data.filter_cols.type.sidebar': extra_params.get('table_data_filter_cols_type_sidebar', table_data_filter_cols_type_sidebar),
            'table.data.filter_cols.tables.value': extra_params.get('table_data_filter_cols_tables_value', table_data_filter_cols_tables_value),
            'table.data.filter_cols.tables.sidebar': extra_params.get('table_data_filter_cols_tables_sidebar', table_data_filter_cols_tables_sidebar),
            'table.data.filter_cols.cname.value': extra_params.get('table_data_filter_cols_cname_value', table_data_filter_cols_cname_value),
            'table.data.filter_cols.cname.sidebar': extra_params.get('table_data_filter_cols_cname_sidebar', table_data_filter_cols_cname_sidebar),
            'table.data.filter_cols.expr.value': extra_params.get('table_data_filter_cols_expr_value', table_data_filter_cols_expr_value),
            'table.data.filter_cols.expr.sidebar': extra_params.get('table_data_filter_cols_expr_sidebar', table_data_filter_cols_expr_sidebar),
            'table.data.filter_cols.value_enum.value': extra_params.get('table_data_filter_cols_value_enum_value', table_data_filter_cols_value_enum_value),
            'table.data.filter_cols.value_enum.sidebar': extra_params.get('table_data_filter_cols_value_enum_sidebar', table_data_filter_cols_value_enum_sidebar),
            'table.data.filter_cols.value_str.value': extra_params.get('table_data_filter_cols_value_str_value', table_data_filter_cols_value_str_value),
            'table.data.filter_cols.value_str.sidebar': extra_params.get('table_data_filter_cols_value_str_sidebar', table_data_filter_cols_value_str_sidebar),
            'table.data.filter_rows.sidebar': extra_params.get('table_data_filter_rows_sidebar', table_data_filter_rows_sidebar),
            'table.data.filter_rows.type.value': extra_params.get('table_data_filter_rows_type_value', table_data_filter_rows_type_value),
            'table.data.filter_rows.type.sidebar': extra_params.get('table_data_filter_rows_type_sidebar', table_data_filter_rows_type_sidebar),
            'table.data.filter_rows.tables.value': extra_params.get('table_data_filter_rows_tables_value', table_data_filter_rows_tables_value),
            'table.data.filter_rows.tables.sidebar': extra_params.get('table_data_filter_rows_tables_sidebar', table_data_filter_rows_tables_sidebar),
            'table.data.filter_rows.cname.value': extra_params.get('table_data_filter_rows_cname_value', table_data_filter_rows_cname_value),
            'table.data.filter_rows.cname.sidebar': extra_params.get('table_data_filter_rows_cname_sidebar', table_data_filter_rows_cname_sidebar),
            'table.data.filter_rows.expr.value': extra_params.get('table_data_filter_rows_expr_value', table_data_filter_rows_expr_value),
            'table.data.filter_rows.expr.sidebar': extra_params.get('table_data_filter_rows_expr_sidebar', table_data_filter_rows_expr_sidebar),
            'table.data.filter_rows.value_enum.value': extra_params.get('table_data_filter_rows_value_enum_value', table_data_filter_rows_value_enum_value),
            'table.data.filter_rows.value_enum.sidebar': extra_params.get('table_data_filter_rows_value_enum_sidebar', table_data_filter_rows_value_enum_sidebar),
            'table.data.filter_rows.value_str.value': extra_params.get('table_data_filter_rows_value_str_value', table_data_filter_rows_value_str_value),
            'table.data.filter_rows.value_str.sidebar': extra_params.get('table_data_filter_rows_value_str_sidebar', table_data_filter_rows_value_str_sidebar),
            'table.category.enabled.value': extra_params.get('table_category_enabled_value', table_category_enabled_value),
            'table.category.enabled.sidebar': extra_params.get('table_category_enabled_sidebar', table_category_enabled_sidebar),
            'table.category.sidebar': extra_params.get('table_category_sidebar', table_category_sidebar),
            'table.category.axis.value': extra_params.get('table_category_axis_value', table_category_axis_value),
            'table.category.axis.sidebar': extra_params.get('table_category_axis_sidebar', table_category_axis_sidebar),
            'table.category.transforms.value': extra_params.get('table_category_transforms_value', table_category_transforms_value),
            'table.category.transforms.sidebar': extra_params.get('table_category_transforms_sidebar', table_category_transforms_sidebar),
            'table.category.columns.category.sidebar': extra_params.get('table_category_columns_category_sidebar', table_category_columns_category_sidebar),
            'table.category.columns.category.table.value': extra_params.get('table_category_columns_category_table_value', table_category_columns_category_table_value),
            'table.category.columns.category.table.sidebar': extra_params.get('table_category_columns_category_table_sidebar', table_category_columns_category_table_sidebar),
            'table.category.columns.category.cname.value': extra_params.get('table_category_columns_category_cname_value', table_category_columns_category_cname_value),
            'table.category.columns.category.cname.sidebar': extra_params.get('table_category_columns_category_cname_sidebar', table_category_columns_category_cname_sidebar),
            'table.category.columns.category.label.value': extra_params.get('table_category_columns_category_label_value', table_category_columns_category_label_value),
            'table.category.columns.category.label.sidebar': extra_params.get('table_category_columns_category_label_sidebar', table_category_columns_category_label_sidebar),
            'table.category.columns.category.scale.value': extra_params.get('table_category_columns_category_scale_value', table_category_columns_category_scale_value),
            'table.category.columns.category.scale.sidebar': extra_params.get('table_category_columns_category_scale_sidebar', table_category_columns_category_scale_sidebar),
            'table.category.columns.category.colorscale': extra_params.get('table_category_columns_category_colorscale', table_category_columns_category_colorscale),
            'table.category.columns.category.is_categorical.value': extra_params.get('table_category_columns_category_is_categorical_value', table_category_columns_category_is_categorical_value),
            'table.category.columns.category.is_categorical.sidebar': extra_params.get('table_category_columns_category_is_categorical_sidebar', table_category_columns_category_is_categorical_sidebar),
            'table.category.filter_rows.sidebar': extra_params.get('table_category_filter_rows_sidebar', table_category_filter_rows_sidebar),
            'table.category.filter_rows.type.value': extra_params.get('table_category_filter_rows_type_value', table_category_filter_rows_type_value),
            'table.category.filter_rows.type.sidebar': extra_params.get('table_category_filter_rows_type_sidebar', table_category_filter_rows_type_sidebar),
            'table.category.filter_rows.tables.value': extra_params.get('table_category_filter_rows_tables_value', table_category_filter_rows_tables_value),
            'table.category.filter_rows.tables.sidebar': extra_params.get('table_category_filter_rows_tables_sidebar', table_category_filter_rows_tables_sidebar),
            'table.category.filter_rows.cname.value': extra_params.get('table_category_filter_rows_cname_value', table_category_filter_rows_cname_value),
            'table.category.filter_rows.cname.sidebar': extra_params.get('table_category_filter_rows_cname_sidebar', table_category_filter_rows_cname_sidebar),
            'table.category.filter_rows.expr.value': extra_params.get('table_category_filter_rows_expr_value', table_category_filter_rows_expr_value),
            'table.category.filter_rows.expr.sidebar': extra_params.get('table_category_filter_rows_expr_sidebar', table_category_filter_rows_expr_sidebar),
            'table.category.filter_rows.value_enum.value': extra_params.get('table_category_filter_rows_value_enum_value', table_category_filter_rows_value_enum_value),
            'table.category.filter_rows.value_enum.sidebar': extra_params.get('table_category_filter_rows_value_enum_sidebar', table_category_filter_rows_value_enum_sidebar),
            'table.category.filter_rows.value_str.value': extra_params.get('table_category_filter_rows_value_str_value', table_category_filter_rows_value_str_value),
            'table.category.filter_rows.value_str.sidebar': extra_params.get('table_category_filter_rows_value_str_sidebar', table_category_filter_rows_value_str_sidebar),
            'variable_options.axis.value': extra_params.get('variable_options_axis_value', variable_options_axis_value),
            'variable_options.axis.sidebar': extra_params.get('variable_options_axis_sidebar', variable_options_axis_sidebar),
            'variable_options.log_values.value': extra_params.get('variable_options_log_values_value', variable_options_log_values_value),
            'variable_options.log_values.sidebar': extra_params.get('variable_options_log_values_sidebar', variable_options_log_values_sidebar),
            'variable_options.sort_by.value': extra_params.get('variable_options_sort_by_value', variable_options_sort_by_value),
            'variable_options.sort_by.sidebar': extra_params.get('variable_options_sort_by_sidebar', variable_options_sort_by_sidebar),
            'category_options.axis.value': extra_params.get('category_options_axis_value', category_options_axis_value),
            'category_options.axis.sidebar': extra_params.get('category_options_axis_sidebar', category_options_axis_sidebar),
            'category_options.sort_by.value': extra_params.get('category_options_sort_by_value', category_options_sort_by_value),
            'category_options.sort_by.sidebar': extra_params.get('category_options_sort_by_sidebar', category_options_sort_by_sidebar),
            'display_options.ncols.value': extra_params.get('display_options_ncols_value', display_options_ncols_value),
            'display_options.ncols.sidebar': extra_params.get('display_options_ncols_sidebar', display_options_ncols_sidebar),
            'display_options.outliers.enabled.value': extra_params.get('display_options_outliers_enabled_value', display_options_outliers_enabled_value),
            'display_options.outliers.enabled.sidebar': extra_params.get('display_options_outliers_enabled_sidebar', display_options_outliers_enabled_sidebar),
            'display_options.title.value': extra_params.get('display_options_title_value', display_options_title_value),
            'display_options.title.sidebar': extra_params.get('display_options_title_sidebar', display_options_title_sidebar),
            'display_options.var_label.value': extra_params.get('display_options_var_label_value', display_options_var_label_value),
            'display_options.var_label.sidebar': extra_params.get('display_options_var_label_sidebar', display_options_var_label_sidebar),
            'display_options.val_label.value': extra_params.get('display_options_val_label_value', display_options_val_label_value),
            'display_options.val_label.sidebar': extra_params.get('display_options_val_label_sidebar', display_options_val_label_sidebar),
            'display_options.height.value': extra_params.get('display_options_height_value', display_options_height_value),
            'display_options.height.sidebar': extra_params.get('display_options_height_sidebar', display_options_height_sidebar),
            'display_options.legend.value': extra_params.get('display_options_legend_value', display_options_legend_value),
            'display_options.legend.sidebar': extra_params.get('display_options_legend_sidebar', display_options_legend_sidebar)
        }
    )
