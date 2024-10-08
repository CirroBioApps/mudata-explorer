# THIS FILE IS AUTOGENERATED

from mudata_explorer.sdk.helpers import collapse_params
from mudata_explorer.helpers.assets import make_process
from mudata import MuData


def summary_stats(
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
    outputs_dest_key_value='summary_stats',
    outputs_dest_key_sidebar=False,
    **extra_params
):
    """
    
    Calculate a variety of summary statistics for the selected data.

    Note that only numerical data may be summarized in this way.

    - count: Number of non-null values
    - prop_valid: Proportion of non-null values
    - prop_valid_rank: Rank ordering of prop_valid (highest first)
    - nunique: Number of unique values
    - nunique_rank: Rank ordering of number of unique values (highest first)
    - median: Median value
    - median_rank: Rank ordering of median value (highest first)
    - mean: Mean value
    - mean_rank: Rank ordering of mean value (highest first)
    - std: Standard deviation
    - min: Minimum value
    - max: Maximum value
    - 25%: 25th percentile
    - 50%: 50th percentile
    - 75%: 75th percentile
    - prop_positive: Proportion of samples with positive values
    - prop_negative: Proportion of samples with negative values

    
    """

    assert isinstance(mdata, MuData), "mdata must be a MuData object"
    extra_params = collapse_params(extra_params)

    # Instantiate the process using all of the parameters
    process = make_process(
        'summary-stats',
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

