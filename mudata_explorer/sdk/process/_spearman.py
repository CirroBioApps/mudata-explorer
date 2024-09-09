# THIS FILE IS AUTOGENERATED

from mudata_explorer.sdk.helpers import collapse_params
from mudata_explorer.helpers.assets import make_process
from mudata import MuData


def spearman(
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
    table_comparitor_sidebar=False,
    table_comparitor_axis_value=None,
    table_comparitor_axis_sidebar=False,
    table_comparitor_transforms_value=[],
    table_comparitor_transforms_sidebar=False,
    table_comparitor_columns_comparitor_sidebar=False,
    table_comparitor_columns_comparitor_table_value=None,
    table_comparitor_columns_comparitor_table_sidebar=False,
    table_comparitor_columns_comparitor_cname_value=None,
    table_comparitor_columns_comparitor_cname_sidebar=False,
    table_comparitor_columns_comparitor_label_value='Comparitor',
    table_comparitor_columns_comparitor_label_sidebar=False,
    table_comparitor_columns_comparitor_scale_value=None,
    table_comparitor_columns_comparitor_scale_sidebar=False,
    table_comparitor_columns_comparitor_colorscale=False,
    table_comparitor_columns_comparitor_is_categorical_value=False,
    table_comparitor_columns_comparitor_is_categorical_sidebar=False,
    table_comparitor_filter_rows_sidebar=False,
    table_comparitor_filter_rows_type_value=None,
    table_comparitor_filter_rows_type_sidebar=False,
    table_comparitor_filter_rows_tables_value=[],
    table_comparitor_filter_rows_tables_sidebar=False,
    table_comparitor_filter_rows_cname_value=None,
    table_comparitor_filter_rows_cname_sidebar=False,
    table_comparitor_filter_rows_expr_value=None,
    table_comparitor_filter_rows_expr_sidebar=False,
    table_comparitor_filter_rows_value_enum_value=None,
    table_comparitor_filter_rows_value_enum_sidebar=False,
    table_comparitor_filter_rows_value_str_value=None,
    table_comparitor_filter_rows_value_str_sidebar=False,
    outputs_dest_key_value='spearman',
    outputs_dest_key_sidebar=False,
    **extra_params
):
    """
    
The Spearman rank-order correlation coefficient is a nonparametric measure of
the strength and direction of association between two ranked variables.
It assesses how well the relationship between two variables can be described
using a monotonic function.
A simple way of describing the analysis is that it is a Pearson correlation
on the ranks of the data (and so it does not assume a linear relationship).

Documentation:
- Wikipedia: [Spearman's rank correlation coefficient](https://en.wikipedia.org/wiki/Spearman%27s_rank_correlation_coefficient)
- Scipy: [scipy.stats.spearmanr](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.spearmanr.html)
    
    """

    assert isinstance(mdata, MuData), "mdata must be a MuData object"
    extra_params = collapse_params(extra_params)

    # Instantiate the process using all of the parameters
    process = make_process(
        'spearman',
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
            'table.comparitor.sidebar': extra_params.get('table_comparitor_sidebar', table_comparitor_sidebar),
            'table.comparitor.axis.value': extra_params.get('table_comparitor_axis_value', table_comparitor_axis_value),
            'table.comparitor.axis.sidebar': extra_params.get('table_comparitor_axis_sidebar', table_comparitor_axis_sidebar),
            'table.comparitor.transforms.value': extra_params.get('table_comparitor_transforms_value', table_comparitor_transforms_value),
            'table.comparitor.transforms.sidebar': extra_params.get('table_comparitor_transforms_sidebar', table_comparitor_transforms_sidebar),
            'table.comparitor.columns.comparitor.sidebar': extra_params.get('table_comparitor_columns_comparitor_sidebar', table_comparitor_columns_comparitor_sidebar),
            'table.comparitor.columns.comparitor.table.value': extra_params.get('table_comparitor_columns_comparitor_table_value', table_comparitor_columns_comparitor_table_value),
            'table.comparitor.columns.comparitor.table.sidebar': extra_params.get('table_comparitor_columns_comparitor_table_sidebar', table_comparitor_columns_comparitor_table_sidebar),
            'table.comparitor.columns.comparitor.cname.value': extra_params.get('table_comparitor_columns_comparitor_cname_value', table_comparitor_columns_comparitor_cname_value),
            'table.comparitor.columns.comparitor.cname.sidebar': extra_params.get('table_comparitor_columns_comparitor_cname_sidebar', table_comparitor_columns_comparitor_cname_sidebar),
            'table.comparitor.columns.comparitor.label.value': extra_params.get('table_comparitor_columns_comparitor_label_value', table_comparitor_columns_comparitor_label_value),
            'table.comparitor.columns.comparitor.label.sidebar': extra_params.get('table_comparitor_columns_comparitor_label_sidebar', table_comparitor_columns_comparitor_label_sidebar),
            'table.comparitor.columns.comparitor.scale.value': extra_params.get('table_comparitor_columns_comparitor_scale_value', table_comparitor_columns_comparitor_scale_value),
            'table.comparitor.columns.comparitor.scale.sidebar': extra_params.get('table_comparitor_columns_comparitor_scale_sidebar', table_comparitor_columns_comparitor_scale_sidebar),
            'table.comparitor.columns.comparitor.colorscale': extra_params.get('table_comparitor_columns_comparitor_colorscale', table_comparitor_columns_comparitor_colorscale),
            'table.comparitor.columns.comparitor.is_categorical.value': extra_params.get('table_comparitor_columns_comparitor_is_categorical_value', table_comparitor_columns_comparitor_is_categorical_value),
            'table.comparitor.columns.comparitor.is_categorical.sidebar': extra_params.get('table_comparitor_columns_comparitor_is_categorical_sidebar', table_comparitor_columns_comparitor_is_categorical_sidebar),
            'table.comparitor.filter_rows.sidebar': extra_params.get('table_comparitor_filter_rows_sidebar', table_comparitor_filter_rows_sidebar),
            'table.comparitor.filter_rows.type.value': extra_params.get('table_comparitor_filter_rows_type_value', table_comparitor_filter_rows_type_value),
            'table.comparitor.filter_rows.type.sidebar': extra_params.get('table_comparitor_filter_rows_type_sidebar', table_comparitor_filter_rows_type_sidebar),
            'table.comparitor.filter_rows.tables.value': extra_params.get('table_comparitor_filter_rows_tables_value', table_comparitor_filter_rows_tables_value),
            'table.comparitor.filter_rows.tables.sidebar': extra_params.get('table_comparitor_filter_rows_tables_sidebar', table_comparitor_filter_rows_tables_sidebar),
            'table.comparitor.filter_rows.cname.value': extra_params.get('table_comparitor_filter_rows_cname_value', table_comparitor_filter_rows_cname_value),
            'table.comparitor.filter_rows.cname.sidebar': extra_params.get('table_comparitor_filter_rows_cname_sidebar', table_comparitor_filter_rows_cname_sidebar),
            'table.comparitor.filter_rows.expr.value': extra_params.get('table_comparitor_filter_rows_expr_value', table_comparitor_filter_rows_expr_value),
            'table.comparitor.filter_rows.expr.sidebar': extra_params.get('table_comparitor_filter_rows_expr_sidebar', table_comparitor_filter_rows_expr_sidebar),
            'table.comparitor.filter_rows.value_enum.value': extra_params.get('table_comparitor_filter_rows_value_enum_value', table_comparitor_filter_rows_value_enum_value),
            'table.comparitor.filter_rows.value_enum.sidebar': extra_params.get('table_comparitor_filter_rows_value_enum_sidebar', table_comparitor_filter_rows_value_enum_sidebar),
            'table.comparitor.filter_rows.value_str.value': extra_params.get('table_comparitor_filter_rows_value_str_value', table_comparitor_filter_rows_value_str_value),
            'table.comparitor.filter_rows.value_str.sidebar': extra_params.get('table_comparitor_filter_rows_value_str_sidebar', table_comparitor_filter_rows_value_str_sidebar),
            'outputs.dest_key.value': extra_params.get('outputs_dest_key_value', outputs_dest_key_value),
            'outputs.dest_key.sidebar': extra_params.get('outputs_dest_key_sidebar', outputs_dest_key_sidebar)
        },
        mdata=mdata,
        params_editable=False
    )

    assert process.params_editable is False, "params_editable must be False"
    assert isinstance(process.mdata, MuData), type(process.mdata)

    # Populate the params for the process
    process.populate_params()

    # Run the process
    process.execute()

