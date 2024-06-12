# THIS FILE IS AUTOGENERATED

from mudata_explorer import app
from muon import MuData
from mudata_explorer.helpers import make_process


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
    outputs_dest_key='summary_stats'
):
    """
    
    Calculate a variety of summary statistics for the selected data.

    Note that only numerical columns may be summarized in this way.

    - count: Number of non-null values in each column
    - prop_valid: Proportion of non-null values in each column
    - nunique: Number of unique values in each column
    - median: Median value of each column
    - mean: Mean value of each column
    - std: Standard deviation of each column
    - min: Minimum value of each column
    - max: Maximum value of each column
    - 25%: 25th percentile of each column
    - 50%: 50th percentile of each column
    - 75%: 75th percentile of each column

    
    """

    assert isinstance(mdata, MuData), "mdata must be a MuData object"

    # Instantiate the process using all of the parameters
    process = make_process(
        'summary-stats',
        params={
            'table.data.axis': table_data_axis,
            'table.data.tables': table_data_tables,
            'table.data.rows_query.query.type': table_data_rows_query_query_type,
            'table.data.rows_query.query.table': table_data_rows_query_query_table,
            'table.data.rows_query.query.cname': table_data_rows_query_query_cname,
            'table.data.rows_query.query.expr': table_data_rows_query_query_expr,
            'table.data.rows_query.query.value': table_data_rows_query_query_value,
            'table.data.cols_query.query.type': table_data_cols_query_query_type,
            'table.data.cols_query.query.table': table_data_cols_query_query_table,
            'table.data.cols_query.query.cname': table_data_cols_query_query_cname,
            'table.data.cols_query.query.expr': table_data_cols_query_query_expr,
            'table.data.cols_query.query.value': table_data_cols_query_query_value,
            'table.data.transforms': table_data_transforms,
            'outputs.dest_key': outputs_dest_key
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

