# THIS FILE IS AUTOGENERATED

from mudata_explorer.sdk.helpers import collapse_params
from mudata_explorer.helpers.add_view import add_view
from muon import MuData


def seaborn_clustermap(
    mdata: MuData,
    data_sidebar=False,
    data_axis_value=None,
    data_axis_sidebar=False,
    data_transforms_value=[],
    data_transforms_sidebar=False,
    data_tables_value=[],
    data_tables_sidebar=False,
    data_filter_cols_sidebar=False,
    data_filter_cols_type_value=None,
    data_filter_cols_type_sidebar=False,
    data_filter_cols_tables_value=[],
    data_filter_cols_tables_sidebar=False,
    data_filter_cols_cname_value=None,
    data_filter_cols_cname_sidebar=False,
    data_filter_cols_expr_value=None,
    data_filter_cols_expr_sidebar=False,
    data_filter_cols_value_enum_value=None,
    data_filter_cols_value_enum_sidebar=False,
    data_filter_cols_value_str_value=None,
    data_filter_cols_value_str_sidebar=False,
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
    z_score_value='None',
    z_score_sidebar=False,
    **extra_params
):
    """
    Display a heatmap with clustered rows and columns using Seaborn.
    """

    assert isinstance(mdata, MuData), "mdata must be a MuData object"
    extra_params = collapse_params(extra_params)

    add_view(
        'seaborn-clustermap',
        mdata,
        params={
            'data.sidebar': extra_params.get('data_sidebar', data_sidebar),
            'data.axis.value': extra_params.get('data_axis_value', data_axis_value),
            'data.axis.sidebar': extra_params.get('data_axis_sidebar', data_axis_sidebar),
            'data.transforms.value': extra_params.get('data_transforms_value', data_transforms_value),
            'data.transforms.sidebar': extra_params.get('data_transforms_sidebar', data_transforms_sidebar),
            'data.tables.value': extra_params.get('data_tables_value', data_tables_value),
            'data.tables.sidebar': extra_params.get('data_tables_sidebar', data_tables_sidebar),
            'data.filter_cols.sidebar': extra_params.get('data_filter_cols_sidebar', data_filter_cols_sidebar),
            'data.filter_cols.type.value': extra_params.get('data_filter_cols_type_value', data_filter_cols_type_value),
            'data.filter_cols.type.sidebar': extra_params.get('data_filter_cols_type_sidebar', data_filter_cols_type_sidebar),
            'data.filter_cols.tables.value': extra_params.get('data_filter_cols_tables_value', data_filter_cols_tables_value),
            'data.filter_cols.tables.sidebar': extra_params.get('data_filter_cols_tables_sidebar', data_filter_cols_tables_sidebar),
            'data.filter_cols.cname.value': extra_params.get('data_filter_cols_cname_value', data_filter_cols_cname_value),
            'data.filter_cols.cname.sidebar': extra_params.get('data_filter_cols_cname_sidebar', data_filter_cols_cname_sidebar),
            'data.filter_cols.expr.value': extra_params.get('data_filter_cols_expr_value', data_filter_cols_expr_value),
            'data.filter_cols.expr.sidebar': extra_params.get('data_filter_cols_expr_sidebar', data_filter_cols_expr_sidebar),
            'data.filter_cols.value_enum.value': extra_params.get('data_filter_cols_value_enum_value', data_filter_cols_value_enum_value),
            'data.filter_cols.value_enum.sidebar': extra_params.get('data_filter_cols_value_enum_sidebar', data_filter_cols_value_enum_sidebar),
            'data.filter_cols.value_str.value': extra_params.get('data_filter_cols_value_str_value', data_filter_cols_value_str_value),
            'data.filter_cols.value_str.sidebar': extra_params.get('data_filter_cols_value_str_sidebar', data_filter_cols_value_str_sidebar),
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
            'z_score.value': extra_params.get('z_score_value', z_score_value),
            'z_score.sidebar': extra_params.get('z_score_sidebar', z_score_sidebar)
        }
    )
