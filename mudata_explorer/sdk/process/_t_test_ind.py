# THIS FILE IS AUTOGENERATED

from mudata_explorer.sdk.helpers import collapse_params
from mudata_explorer.helpers.assets import make_process
from muon import MuData


def t_test_ind(
    mdata: MuData,
    data_axis=0,
    data_tables=[],
    data_rows_query_query_type='',
    data_rows_query_query_table='',
    data_rows_query_query_cname='',
    data_rows_query_query_expr='',
    data_rows_query_query_value='',
    data_cols_query_query_type='',
    data_cols_query_query_table='',
    data_cols_query_query_cname='',
    data_cols_query_query_expr='',
    data_cols_query_query_value='',
    data_transforms=[],
    group_axis=0,
    group_group_table=None,
    group_group_cname=None,
    group_group_label='Group',
    group_rows_query_query_type='',
    group_rows_query_query_table='',
    group_rows_query_query_cname='',
    group_rows_query_query_expr='',
    group_rows_query_query_value='',
    group_cols_query_query_type='',
    group_cols_query_query_table='',
    group_cols_query_query_cname='',
    group_cols_query_query_expr='',
    group_cols_query_query_value='',
    group_transforms=[],
    dest_key='t_test_ind',
    **extra_params
):
    """
    
    The independent samples t-test is used to determine if there is a significant
    difference between the means of two independent groups. It assumes that the
    data is normally distributed and that the variances of the two groups are equal.

    > Make sure that the 'Group' column only contains two unique values.

    Results will include:

    - pvalue: The p-value of the t-test
    - neg_log10_pvalue: The negative log10 of the p-value
    - t_statistic: The t-statistic of the t-test
    - median_diff: The difference in medians between the two groups

    
    """

    assert isinstance(mdata, MuData), "mdata must be a MuData object"
    extra_params = collapse_params(extra_params)

    # Instantiate the process using all of the parameters
    process = make_process(
        't-test-ind',
        params={
            'data.axis': extra_params.get('data_axis', data_axis),
            'data.tables': extra_params.get('data_tables', data_tables),
            'data.rows_query.query.type': extra_params.get('data_rows_query_query_type', data_rows_query_query_type),
            'data.rows_query.query.table': extra_params.get('data_rows_query_query_table', data_rows_query_query_table),
            'data.rows_query.query.cname': extra_params.get('data_rows_query_query_cname', data_rows_query_query_cname),
            'data.rows_query.query.expr': extra_params.get('data_rows_query_query_expr', data_rows_query_query_expr),
            'data.rows_query.query.value': extra_params.get('data_rows_query_query_value', data_rows_query_query_value),
            'data.cols_query.query.type': extra_params.get('data_cols_query_query_type', data_cols_query_query_type),
            'data.cols_query.query.table': extra_params.get('data_cols_query_query_table', data_cols_query_query_table),
            'data.cols_query.query.cname': extra_params.get('data_cols_query_query_cname', data_cols_query_query_cname),
            'data.cols_query.query.expr': extra_params.get('data_cols_query_query_expr', data_cols_query_query_expr),
            'data.cols_query.query.value': extra_params.get('data_cols_query_query_value', data_cols_query_query_value),
            'data.transforms': extra_params.get('data_transforms', data_transforms),
            'group.axis': extra_params.get('group_axis', group_axis),
            'group.group.table': extra_params.get('group_group_table', group_group_table),
            'group.group.cname': extra_params.get('group_group_cname', group_group_cname),
            'group.group.label': extra_params.get('group_group_label', group_group_label),
            'group.rows_query.query.type': extra_params.get('group_rows_query_query_type', group_rows_query_query_type),
            'group.rows_query.query.table': extra_params.get('group_rows_query_query_table', group_rows_query_query_table),
            'group.rows_query.query.cname': extra_params.get('group_rows_query_query_cname', group_rows_query_query_cname),
            'group.rows_query.query.expr': extra_params.get('group_rows_query_query_expr', group_rows_query_query_expr),
            'group.rows_query.query.value': extra_params.get('group_rows_query_query_value', group_rows_query_query_value),
            'group.cols_query.query.type': extra_params.get('group_cols_query_query_type', group_cols_query_query_type),
            'group.cols_query.query.table': extra_params.get('group_cols_query_query_table', group_cols_query_query_table),
            'group.cols_query.query.cname': extra_params.get('group_cols_query_query_cname', group_cols_query_query_cname),
            'group.cols_query.query.expr': extra_params.get('group_cols_query_query_expr', group_cols_query_query_expr),
            'group.cols_query.query.value': extra_params.get('group_cols_query_query_value', group_cols_query_query_value),
            'group.transforms': extra_params.get('group_transforms', group_transforms),
            'dest_key': extra_params.get('dest_key', dest_key)
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

