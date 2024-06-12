# THIS FILE IS AUTOGENERATED

from mudata_explorer import app
from muon import MuData
from mudata_explorer.helpers import make_process


def decision_tree_classifier(
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
    table_predictor_axis=0,
    table_predictor_predictor_table=None,
    table_predictor_predictor_cname=None,
    table_predictor_predictor_label='Predictor',
    table_predictor_rows_query_query_type='',
    table_predictor_rows_query_query_table='',
    table_predictor_rows_query_query_cname='',
    table_predictor_rows_query_query_expr='',
    table_predictor_rows_query_query_value='',
    table_predictor_cols_query_query_type='',
    table_predictor_cols_query_query_table='',
    table_predictor_cols_query_query_cname='',
    table_predictor_cols_query_query_expr='',
    table_predictor_cols_query_query_value='',
    table_predictor_transforms=[],
    model_params_train_prop=0.5,
    model_params_criterion='gini',
    model_params_max_depth=0,
    model_params_min_samples_split=2,
    model_params_random_state=0,
    outputs_dest_key='decision_tree_classifier'
):
    """
    
    Predict the value of a target variable by learning
    simple decision rules inferred from the data features.
    
    This analysis is used to predict _categorical_ variables.
    
    The Decision Tree Classifier is implemented using the
    scikit-learn library, and more complete information can
    be found on the [scikit-learn documentation](https://scikit-learn.org/stable/modules/generated/sklearn.tree.DecisionTreeClassifier.html).

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

    # Instantiate the process using all of the parameters
    process = make_process(
        'decision_tree_classifier',
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
            'table.predictor.axis': table_predictor_axis,
            'table.predictor.predictor.table': table_predictor_predictor_table,
            'table.predictor.predictor.cname': table_predictor_predictor_cname,
            'table.predictor.predictor.label': table_predictor_predictor_label,
            'table.predictor.rows_query.query.type': table_predictor_rows_query_query_type,
            'table.predictor.rows_query.query.table': table_predictor_rows_query_query_table,
            'table.predictor.rows_query.query.cname': table_predictor_rows_query_query_cname,
            'table.predictor.rows_query.query.expr': table_predictor_rows_query_query_expr,
            'table.predictor.rows_query.query.value': table_predictor_rows_query_query_value,
            'table.predictor.cols_query.query.type': table_predictor_cols_query_query_type,
            'table.predictor.cols_query.query.table': table_predictor_cols_query_query_table,
            'table.predictor.cols_query.query.cname': table_predictor_cols_query_query_cname,
            'table.predictor.cols_query.query.expr': table_predictor_cols_query_query_expr,
            'table.predictor.cols_query.query.value': table_predictor_cols_query_query_value,
            'table.predictor.transforms': table_predictor_transforms,
            'model_params.train_prop': model_params_train_prop,
            'model_params.criterion': model_params_criterion,
            'model_params.max_depth': model_params_max_depth,
            'model_params.min_samples_split': model_params_min_samples_split,
            'model_params.random_state': model_params_random_state,
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

