# THIS FILE IS AUTOGENERATED

from mudata_explorer.sdk.helpers import collapse_params
from mudata_explorer.helpers.add_view import add_view
from muon import MuData


def plotly_stacked_bars(
    mdata: MuData,
    table_data_sidebar=False,
    table_data_axis_value=0,
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
    table_category_axis_value=0,
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
    formatting_max_features_enabled_value=True,
    formatting_max_features_enabled_sidebar=False,
    formatting_max_features_value=0,
    formatting_max_features_sidebar=False,
    formatting_sort_cols_by_value='Labels',
    formatting_sort_cols_by_sidebar=True,
    formatting_sort_rows_by_value='Labels',
    formatting_sort_rows_by_sidebar=True,
    formatting_title_value='',
    formatting_title_sidebar=True,
    formatting_yaxis_title_value='',
    formatting_yaxis_title_sidebar=False,
    formatting_feature_label_value='',
    formatting_feature_label_sidebar=False,
    formatting_max_y_value=0,
    formatting_max_y_sidebar=False,
    **extra_params
):
    """
    
Show the values in a collection of columns, optionally
annotated by a single column which contains categories.
    
    """

    assert isinstance(mdata, MuData), "mdata must be a MuData object"
    extra_params = collapse_params(extra_params)

    add_view(
        'plotly-stacked-bars',
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
            'formatting.max_features.enabled.value': extra_params.get('formatting_max_features_enabled_value', formatting_max_features_enabled_value),
            'formatting.max_features.enabled.sidebar': extra_params.get('formatting_max_features_enabled_sidebar', formatting_max_features_enabled_sidebar),
            'formatting.max_features.value': extra_params.get('formatting_max_features_value', formatting_max_features_value),
            'formatting.max_features.sidebar': extra_params.get('formatting_max_features_sidebar', formatting_max_features_sidebar),
            'formatting.sort_cols_by.value': extra_params.get('formatting_sort_cols_by_value', formatting_sort_cols_by_value),
            'formatting.sort_cols_by.sidebar': extra_params.get('formatting_sort_cols_by_sidebar', formatting_sort_cols_by_sidebar),
            'formatting.sort_rows_by.value': extra_params.get('formatting_sort_rows_by_value', formatting_sort_rows_by_value),
            'formatting.sort_rows_by.sidebar': extra_params.get('formatting_sort_rows_by_sidebar', formatting_sort_rows_by_sidebar),
            'formatting.title.value': extra_params.get('formatting_title_value', formatting_title_value),
            'formatting.title.sidebar': extra_params.get('formatting_title_sidebar', formatting_title_sidebar),
            'formatting.yaxis_title.value': extra_params.get('formatting_yaxis_title_value', formatting_yaxis_title_value),
            'formatting.yaxis_title.sidebar': extra_params.get('formatting_yaxis_title_sidebar', formatting_yaxis_title_sidebar),
            'formatting.feature_label.value': extra_params.get('formatting_feature_label_value', formatting_feature_label_value),
            'formatting.feature_label.sidebar': extra_params.get('formatting_feature_label_sidebar', formatting_feature_label_sidebar),
            'formatting.max_y.value': extra_params.get('formatting_max_y_value', formatting_max_y_value),
            'formatting.max_y.sidebar': extra_params.get('formatting_max_y_sidebar', formatting_max_y_sidebar)
        }
    )
