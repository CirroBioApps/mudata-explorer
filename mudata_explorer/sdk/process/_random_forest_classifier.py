# THIS FILE IS AUTOGENERATED

from mudata_explorer.sdk.helpers import collapse_params
from mudata_explorer.helpers.assets import make_process
from muon import MuData


def random_forest_classifier(
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
    table_predictor_sidebar=False,
    table_predictor_axis_value=None,
    table_predictor_axis_sidebar=False,
    table_predictor_transforms_value=[],
    table_predictor_transforms_sidebar=False,
    table_predictor_columns_predictor_sidebar=False,
    table_predictor_columns_predictor_table_value=None,
    table_predictor_columns_predictor_table_sidebar=False,
    table_predictor_columns_predictor_cname_value=None,
    table_predictor_columns_predictor_cname_sidebar=False,
    table_predictor_columns_predictor_label_value='Predictor',
    table_predictor_columns_predictor_label_sidebar=False,
    table_predictor_columns_predictor_scale_value=None,
    table_predictor_columns_predictor_scale_sidebar=False,
    table_predictor_columns_predictor_colorscale=False,
    table_predictor_columns_predictor_is_categorical_value=False,
    table_predictor_columns_predictor_is_categorical_sidebar=False,
    table_predictor_filter_rows_sidebar=False,
    table_predictor_filter_rows_type_value=None,
    table_predictor_filter_rows_type_sidebar=False,
    table_predictor_filter_rows_tables_value=[],
    table_predictor_filter_rows_tables_sidebar=False,
    table_predictor_filter_rows_cname_value=None,
    table_predictor_filter_rows_cname_sidebar=False,
    table_predictor_filter_rows_expr_value=None,
    table_predictor_filter_rows_expr_sidebar=False,
    table_predictor_filter_rows_value_enum_value=None,
    table_predictor_filter_rows_value_enum_sidebar=False,
    table_predictor_filter_rows_value_str_value=None,
    table_predictor_filter_rows_value_str_sidebar=False,
    model_params_train_prop_value=0.5,
    model_params_train_prop_sidebar=False,
    model_params_balance_groups_value='Use All Data',
    model_params_balance_groups_sidebar=False,
    model_params_n_estimators_value=100,
    model_params_n_estimators_sidebar=False,
    model_params_criterion_value='gini',
    model_params_criterion_sidebar=False,
    model_params_max_depth_value=0,
    model_params_max_depth_sidebar=False,
    model_params_min_samples_split_value=2,
    model_params_min_samples_split_sidebar=False,
    model_params_random_state_value=0,
    model_params_random_state_sidebar=False,
    outputs_dest_key_value='random_forest_classifier',
    outputs_dest_key_sidebar=False,
    **extra_params
):
    """
    
    Predict the value of a target variable building a large
    number of decision tree classifiers each using
    simple decision rules inferred from the data features.

    This analysis is used to predict _categorical_ variables.
    
    The Random Forest Classifier is implemented using the
    scikit-learn library, and more complete information can
    be found on the [scikit-learn documentation](https://scikit-learn.org/stable/modules/generated/sklearn.ensemble.RandomForestClassifier.html).

    The performance of the model is evaluated by comparing the
    accuracy of classification on the held out testing set
    against the accuracy observed under random permutation
    of the labels in that testing set.

    The outputs of the analysis include:

    - A table of the actual and predicted values of the target variable,
    which includes a column indicating whether the observation was
    part of the training or testing set.
    - A table of the weights of each feature in the model.
    
    """

    assert isinstance(mdata, MuData), "mdata must be a MuData object"
    extra_params = collapse_params(extra_params)

    # Instantiate the process using all of the parameters
    process = make_process(
        'random_forest_classifier',
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
            'table.predictor.sidebar': extra_params.get('table_predictor_sidebar', table_predictor_sidebar),
            'table.predictor.axis.value': extra_params.get('table_predictor_axis_value', table_predictor_axis_value),
            'table.predictor.axis.sidebar': extra_params.get('table_predictor_axis_sidebar', table_predictor_axis_sidebar),
            'table.predictor.transforms.value': extra_params.get('table_predictor_transforms_value', table_predictor_transforms_value),
            'table.predictor.transforms.sidebar': extra_params.get('table_predictor_transforms_sidebar', table_predictor_transforms_sidebar),
            'table.predictor.columns.predictor.sidebar': extra_params.get('table_predictor_columns_predictor_sidebar', table_predictor_columns_predictor_sidebar),
            'table.predictor.columns.predictor.table.value': extra_params.get('table_predictor_columns_predictor_table_value', table_predictor_columns_predictor_table_value),
            'table.predictor.columns.predictor.table.sidebar': extra_params.get('table_predictor_columns_predictor_table_sidebar', table_predictor_columns_predictor_table_sidebar),
            'table.predictor.columns.predictor.cname.value': extra_params.get('table_predictor_columns_predictor_cname_value', table_predictor_columns_predictor_cname_value),
            'table.predictor.columns.predictor.cname.sidebar': extra_params.get('table_predictor_columns_predictor_cname_sidebar', table_predictor_columns_predictor_cname_sidebar),
            'table.predictor.columns.predictor.label.value': extra_params.get('table_predictor_columns_predictor_label_value', table_predictor_columns_predictor_label_value),
            'table.predictor.columns.predictor.label.sidebar': extra_params.get('table_predictor_columns_predictor_label_sidebar', table_predictor_columns_predictor_label_sidebar),
            'table.predictor.columns.predictor.scale.value': extra_params.get('table_predictor_columns_predictor_scale_value', table_predictor_columns_predictor_scale_value),
            'table.predictor.columns.predictor.scale.sidebar': extra_params.get('table_predictor_columns_predictor_scale_sidebar', table_predictor_columns_predictor_scale_sidebar),
            'table.predictor.columns.predictor.colorscale': extra_params.get('table_predictor_columns_predictor_colorscale', table_predictor_columns_predictor_colorscale),
            'table.predictor.columns.predictor.is_categorical.value': extra_params.get('table_predictor_columns_predictor_is_categorical_value', table_predictor_columns_predictor_is_categorical_value),
            'table.predictor.columns.predictor.is_categorical.sidebar': extra_params.get('table_predictor_columns_predictor_is_categorical_sidebar', table_predictor_columns_predictor_is_categorical_sidebar),
            'table.predictor.filter_rows.sidebar': extra_params.get('table_predictor_filter_rows_sidebar', table_predictor_filter_rows_sidebar),
            'table.predictor.filter_rows.type.value': extra_params.get('table_predictor_filter_rows_type_value', table_predictor_filter_rows_type_value),
            'table.predictor.filter_rows.type.sidebar': extra_params.get('table_predictor_filter_rows_type_sidebar', table_predictor_filter_rows_type_sidebar),
            'table.predictor.filter_rows.tables.value': extra_params.get('table_predictor_filter_rows_tables_value', table_predictor_filter_rows_tables_value),
            'table.predictor.filter_rows.tables.sidebar': extra_params.get('table_predictor_filter_rows_tables_sidebar', table_predictor_filter_rows_tables_sidebar),
            'table.predictor.filter_rows.cname.value': extra_params.get('table_predictor_filter_rows_cname_value', table_predictor_filter_rows_cname_value),
            'table.predictor.filter_rows.cname.sidebar': extra_params.get('table_predictor_filter_rows_cname_sidebar', table_predictor_filter_rows_cname_sidebar),
            'table.predictor.filter_rows.expr.value': extra_params.get('table_predictor_filter_rows_expr_value', table_predictor_filter_rows_expr_value),
            'table.predictor.filter_rows.expr.sidebar': extra_params.get('table_predictor_filter_rows_expr_sidebar', table_predictor_filter_rows_expr_sidebar),
            'table.predictor.filter_rows.value_enum.value': extra_params.get('table_predictor_filter_rows_value_enum_value', table_predictor_filter_rows_value_enum_value),
            'table.predictor.filter_rows.value_enum.sidebar': extra_params.get('table_predictor_filter_rows_value_enum_sidebar', table_predictor_filter_rows_value_enum_sidebar),
            'table.predictor.filter_rows.value_str.value': extra_params.get('table_predictor_filter_rows_value_str_value', table_predictor_filter_rows_value_str_value),
            'table.predictor.filter_rows.value_str.sidebar': extra_params.get('table_predictor_filter_rows_value_str_sidebar', table_predictor_filter_rows_value_str_sidebar),
            'model_params.train_prop.value': extra_params.get('model_params_train_prop_value', model_params_train_prop_value),
            'model_params.train_prop.sidebar': extra_params.get('model_params_train_prop_sidebar', model_params_train_prop_sidebar),
            'model_params.balance_groups.value': extra_params.get('model_params_balance_groups_value', model_params_balance_groups_value),
            'model_params.balance_groups.sidebar': extra_params.get('model_params_balance_groups_sidebar', model_params_balance_groups_sidebar),
            'model_params.n_estimators.value': extra_params.get('model_params_n_estimators_value', model_params_n_estimators_value),
            'model_params.n_estimators.sidebar': extra_params.get('model_params_n_estimators_sidebar', model_params_n_estimators_sidebar),
            'model_params.criterion.value': extra_params.get('model_params_criterion_value', model_params_criterion_value),
            'model_params.criterion.sidebar': extra_params.get('model_params_criterion_sidebar', model_params_criterion_sidebar),
            'model_params.max_depth.value': extra_params.get('model_params_max_depth_value', model_params_max_depth_value),
            'model_params.max_depth.sidebar': extra_params.get('model_params_max_depth_sidebar', model_params_max_depth_sidebar),
            'model_params.min_samples_split.value': extra_params.get('model_params_min_samples_split_value', model_params_min_samples_split_value),
            'model_params.min_samples_split.sidebar': extra_params.get('model_params_min_samples_split_sidebar', model_params_min_samples_split_sidebar),
            'model_params.random_state.value': extra_params.get('model_params_random_state_value', model_params_random_state_value),
            'model_params.random_state.sidebar': extra_params.get('model_params_random_state_sidebar', model_params_random_state_sidebar),
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

