# THIS FILE IS AUTOGENERATED

from mudata_explorer import app
from mudata_explorer.sdk.helpers import collapse_params
from mudata_explorer.helpers.assets import make_process
from muon import MuData


def summary_stats(
    mdata: MuData,
    table_data_axis=0,
    table_data_tables=[],
    table_data_rows_query_query_type='',
    table_data_rows_query_query_table='',
    table_data_rows_query_query_cname='',
    table_data_rows_query_query_expr='',
    table_data_rows_query_query_value='',
    table_data_cols_query_query_type='',
    table_data_cols_query_query_table='',
    table_data_cols_query_query_cname='',
    table_data_cols_query_query_expr='',
    table_data_cols_query_query_value='',
    table_data_transforms=[],
    outputs_dest_key='summary_stats',
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
            'table.data.axis': extra_params.get('table_data_axis', table_data_axis),
            'table.data.tables': extra_params.get('table_data_tables', table_data_tables),
            'table.data.rows_query.query.type': extra_params.get('table_data_rows_query_query_type', table_data_rows_query_query_type),
            'table.data.rows_query.query.table': extra_params.get('table_data_rows_query_query_table', table_data_rows_query_query_table),
            'table.data.rows_query.query.cname': extra_params.get('table_data_rows_query_query_cname', table_data_rows_query_query_cname),
            'table.data.rows_query.query.expr': extra_params.get('table_data_rows_query_query_expr', table_data_rows_query_query_expr),
            'table.data.rows_query.query.value': extra_params.get('table_data_rows_query_query_value', table_data_rows_query_query_value),
            'table.data.cols_query.query.type': extra_params.get('table_data_cols_query_query_type', table_data_cols_query_query_type),
            'table.data.cols_query.query.table': extra_params.get('table_data_cols_query_query_table', table_data_cols_query_query_table),
            'table.data.cols_query.query.cname': extra_params.get('table_data_cols_query_query_cname', table_data_cols_query_query_cname),
            'table.data.cols_query.query.expr': extra_params.get('table_data_cols_query_query_expr', table_data_cols_query_query_expr),
            'table.data.cols_query.query.value': extra_params.get('table_data_cols_query_query_value', table_data_cols_query_query_value),
            'table.data.transforms': extra_params.get('table_data_transforms', table_data_transforms),
            'outputs.dest_key': extra_params.get('outputs_dest_key', outputs_dest_key)
        },
        mdata=mdata,
        params_editable=False
    )

    assert process.params_editable is False, "params_editable must be False"
    assert isinstance(process.mdata, MuData), type(process.mdata)

    # Get the data from the object
    process.get_data()

    # Run the process
    process.execute()

