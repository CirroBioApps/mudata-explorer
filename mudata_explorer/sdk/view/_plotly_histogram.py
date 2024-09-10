# THIS FILE IS AUTOGENERATED

from mudata_explorer.sdk.helpers import collapse_params
from mudata_explorer.helpers.add_view import add_view
from mudata import MuData


def plotly_histogram(
    mdata: MuData,
    data_sidebar=False,
    data_axis_value=None,
    data_axis_sidebar=False,
    data_transforms_value=[],
    data_transforms_sidebar=False,
    data_columns_value_sidebar=False,
    data_columns_value_table_value=None,
    data_columns_value_table_sidebar=False,
    data_columns_value_cname_value=None,
    data_columns_value_cname_sidebar=False,
    data_columns_value_label_value='Value',
    data_columns_value_label_sidebar=False,
    data_columns_value_scale_value=None,
    data_columns_value_scale_sidebar=False,
    data_columns_value_colorscale=False,
    data_columns_value_is_categorical_value=False,
    data_columns_value_is_categorical_sidebar=False,
    data_columns_grouping_enabled_value=True,
    data_columns_grouping_enabled_sidebar=False,
    data_columns_grouping_sidebar=False,
    data_columns_grouping_table_value=None,
    data_columns_grouping_table_sidebar=False,
    data_columns_grouping_cname_value=None,
    data_columns_grouping_cname_sidebar=False,
    data_columns_grouping_label_value='Grouping',
    data_columns_grouping_label_sidebar=False,
    data_columns_grouping_scale_value=None,
    data_columns_grouping_scale_sidebar=False,
    data_columns_grouping_colorscale=True,
    data_columns_grouping_is_categorical_value=True,
    data_columns_grouping_is_categorical_sidebar=False,
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
    scale_options_log_y_value=None,
    scale_options_log_y_sidebar=True,
    formatting_title_value='',
    formatting_title_sidebar=True,
    formatting_legend_value=None,
    formatting_legend_sidebar=False,
    formatting_nbins_value=20,
    formatting_nbins_sidebar=True,
    statistics_compare_groups_value='Disabled',
    statistics_compare_groups_sidebar=True,
    **extra_params
):
    """
    
Display a distribution of values as a frequency histogram using Plotly.

Optionally include a grouping column to display multiple histograms
which are overlaid on the same plot using different colors.
    
    """

    assert isinstance(mdata, MuData), "mdata must be a MuData object"
    extra_params = collapse_params(extra_params)

    add_view(
        'plotly-histogram',
        mdata,
        params={
            'data.sidebar': extra_params.get('data_sidebar', data_sidebar),
            'data.axis.value': extra_params.get('data_axis_value', data_axis_value),
            'data.axis.sidebar': extra_params.get('data_axis_sidebar', data_axis_sidebar),
            'data.transforms.value': extra_params.get('data_transforms_value', data_transforms_value),
            'data.transforms.sidebar': extra_params.get('data_transforms_sidebar', data_transforms_sidebar),
            'data.columns.value.sidebar': extra_params.get('data_columns_value_sidebar', data_columns_value_sidebar),
            'data.columns.value.table.value': extra_params.get('data_columns_value_table_value', data_columns_value_table_value),
            'data.columns.value.table.sidebar': extra_params.get('data_columns_value_table_sidebar', data_columns_value_table_sidebar),
            'data.columns.value.cname.value': extra_params.get('data_columns_value_cname_value', data_columns_value_cname_value),
            'data.columns.value.cname.sidebar': extra_params.get('data_columns_value_cname_sidebar', data_columns_value_cname_sidebar),
            'data.columns.value.label.value': extra_params.get('data_columns_value_label_value', data_columns_value_label_value),
            'data.columns.value.label.sidebar': extra_params.get('data_columns_value_label_sidebar', data_columns_value_label_sidebar),
            'data.columns.value.scale.value': extra_params.get('data_columns_value_scale_value', data_columns_value_scale_value),
            'data.columns.value.scale.sidebar': extra_params.get('data_columns_value_scale_sidebar', data_columns_value_scale_sidebar),
            'data.columns.value.colorscale': extra_params.get('data_columns_value_colorscale', data_columns_value_colorscale),
            'data.columns.value.is_categorical.value': extra_params.get('data_columns_value_is_categorical_value', data_columns_value_is_categorical_value),
            'data.columns.value.is_categorical.sidebar': extra_params.get('data_columns_value_is_categorical_sidebar', data_columns_value_is_categorical_sidebar),
            'data.columns.grouping.enabled.value': extra_params.get('data_columns_grouping_enabled_value', data_columns_grouping_enabled_value),
            'data.columns.grouping.enabled.sidebar': extra_params.get('data_columns_grouping_enabled_sidebar', data_columns_grouping_enabled_sidebar),
            'data.columns.grouping.sidebar': extra_params.get('data_columns_grouping_sidebar', data_columns_grouping_sidebar),
            'data.columns.grouping.table.value': extra_params.get('data_columns_grouping_table_value', data_columns_grouping_table_value),
            'data.columns.grouping.table.sidebar': extra_params.get('data_columns_grouping_table_sidebar', data_columns_grouping_table_sidebar),
            'data.columns.grouping.cname.value': extra_params.get('data_columns_grouping_cname_value', data_columns_grouping_cname_value),
            'data.columns.grouping.cname.sidebar': extra_params.get('data_columns_grouping_cname_sidebar', data_columns_grouping_cname_sidebar),
            'data.columns.grouping.label.value': extra_params.get('data_columns_grouping_label_value', data_columns_grouping_label_value),
            'data.columns.grouping.label.sidebar': extra_params.get('data_columns_grouping_label_sidebar', data_columns_grouping_label_sidebar),
            'data.columns.grouping.scale.value': extra_params.get('data_columns_grouping_scale_value', data_columns_grouping_scale_value),
            'data.columns.grouping.scale.sidebar': extra_params.get('data_columns_grouping_scale_sidebar', data_columns_grouping_scale_sidebar),
            'data.columns.grouping.colorscale': extra_params.get('data_columns_grouping_colorscale', data_columns_grouping_colorscale),
            'data.columns.grouping.is_categorical.value': extra_params.get('data_columns_grouping_is_categorical_value', data_columns_grouping_is_categorical_value),
            'data.columns.grouping.is_categorical.sidebar': extra_params.get('data_columns_grouping_is_categorical_sidebar', data_columns_grouping_is_categorical_sidebar),
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
            'scale_options.log_y.value': extra_params.get('scale_options_log_y_value', scale_options_log_y_value),
            'scale_options.log_y.sidebar': extra_params.get('scale_options_log_y_sidebar', scale_options_log_y_sidebar),
            'formatting.title.value': extra_params.get('formatting_title_value', formatting_title_value),
            'formatting.title.sidebar': extra_params.get('formatting_title_sidebar', formatting_title_sidebar),
            'formatting.legend.value': extra_params.get('formatting_legend_value', formatting_legend_value),
            'formatting.legend.sidebar': extra_params.get('formatting_legend_sidebar', formatting_legend_sidebar),
            'formatting.nbins.value': extra_params.get('formatting_nbins_value', formatting_nbins_value),
            'formatting.nbins.sidebar': extra_params.get('formatting_nbins_sidebar', formatting_nbins_sidebar),
            'statistics.compare_groups.value': extra_params.get('statistics_compare_groups_value', statistics_compare_groups_value),
            'statistics.compare_groups.sidebar': extra_params.get('statistics_compare_groups_sidebar', statistics_compare_groups_sidebar)
        }
    )
