import pandas as pd
import plotly.express as px
from plotly import io
from sklearn import tree
from mudata_explorer.base.process import Process


class DecisionTreeClassifier(Process):

    type = "decision_tree_classifier"
    name = "Decision Tree Classifier"
    help_text = """
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
    """ # noqa
    category = "Machine Learning"
    schema = {
        "table": {
            "type": "object",
            "label": "Data Table",
            "properties": {
                "data": {
                    "type": "dataframe",
                    "label": "Data",
                    "help": "Select the data used to train and test the model",
                    "select_columns": True,
                    "query": "",
                },
                "predictor": {
                    "label": "Predictor",
                    "help": "Select the column containing values to predict",
                    "type": "dataframe",
                    "columns": {"predictor": {"label": "Predictor"}}
                }
            }
        },
        "model_params": {
            "type": "object",
            "label": "Model Parameters",
            "properties": {
                "train_prop": {
                    "type": "float",
                    "min_value": 0.,
                    "max_value": 1.,
                    "default": 0.5,
                    "label": "Training Proportion",
                    "help": "Proportion of the data to use for training."
                },
                "criterion": {
                    "type": "string",
                    "enum": ["gini", "entropy", "log_loss"],
                    "default": "gini",
                    "label": "Criterion",
                    "help": "Function to measure the quality of a split."
                },
                "max_depth": {
                    "type": "integer",
                    "min_value": 0,
                    "default": 0,
                    "label": "Max Depth",
                    "help": "The maximum depth of the tree (0 for unlimited).",
                },
                "min_samples_split": {
                    "type": "integer",
                    "min_value": 2,
                    "default": 2,
                    "label": "Min Samples Split",
                    "help": """
                    The minimum number of samples required to split
                    an internal node.
                    """
                },
                "random_state": {
                    "type": "integer",
                    "default": 0,
                    "label": "Random State",
                    "help": "Controls the randomness of the estimator."
                }
            }
        },
        "outputs": {
            "type": "object",
            "label": "Outputs",
            "properties": {
                "dest_key": {
                    "type": "string",
                    "default": "decision_tree_classifier",
                    "label": "Label to use for results",
                    "help": """
                    Key to use when saving the output
                    """
                }
            }
        }
    }
    outputs = {
        "assignments": {
            "type": pd.DataFrame,
            "label": "Model Assignments",
            "desc": "Results of the model assigning labels to the data.",
            "modality": "table.data.tables",
            "axis": "table.data.axis",
            "attr": "outputs.dest_key"
        },
        "weights": {
            "type": pd.DataFrame,
            "label": "Feature Weights",
            "desc": "Weights of each feature in the model.",
            "modality": "table.data.tables",
            "axis": "table.data.axis.T",
            "attr": "outputs.dest_key"
        }
    }

    def execute(self):

        df: pd.DataFrame = self.params["table.data.dataframe"]
        df = df.dropna()
        assert df.shape[0] > 0, "Null values in all rows."

        pred: pd.Series = (
            self.params
            ["table.predictor.dataframe"]
            ["predictor"]
            .dropna()
        )

        # Get the set of observations which have both data and predictor
        shared = df.index.intersection(pred.index)

        assert len(shared) > 0, "No shared rows between data and predictor"

        df = df.loc[shared]
        pred = pred.loc[shared]

        train_prop = self.params["model_params.train_prop"]
        criterion = self.params["model_params.criterion"]
        max_depth = self.params["model_params.max_depth"]
        min_samples_split = self.params["model_params.min_samples_split"]
        random_state = self.params["model_params.random_state"]

        # Get the training and testing sets
        train_ix = df.sample(frac=train_prop, random_state=random_state).index
        test_ix = df.index.difference(train_ix)

        assert train_ix.shape[0] > 0, "No training data - increase train_prop"
        assert test_ix.shape[0] > 0, "No testing data - decrease train_prop"

        # Run the DecisionTreeClassifier
        clf = tree.DecisionTreeClassifier(
            criterion=criterion,
            max_depth=(
                max_depth
                if max_depth > 0
                else None
            ),
            min_samples_split=min_samples_split,
            random_state=random_state
        )
        clf = clf.fit(df.loc[train_ix], pred.loc[train_ix])

        # Compute the performance of the model by comparing the
        # accuracy of classification on the held out testing set
        # against the accuracy observed under random permutation
        n_reps = 1000
        null_dist = [
            clf.score(
                df.loc[test_ix],
                pred.loc[test_ix].sample(frac=1, random_state=i)
            )
            for i in range(n_reps)
        ]
        score = clf.score(df.loc[test_ix], pred.loc[test_ix])

        fig = px.histogram(
            x=null_dist,
            labels={"x": "Accuracy"},
            title="Model Performance"
        )
        fig.add_vline(
            x=score,
            line_dash="dash",
            line_color="grey"
        )
        # Annotate the figure with text indicating the score value
        perc = sum([n < score for n in null_dist]) / n_reps
        fig.add_annotation(
            x=score,
            xanchor="left",
            y=-20,
            text=f"Model Score: {score:.2f} (Above {int(100 * perc)}% of null distribution)", # noqa
            showarrow=False
        )
        figures = [io.to_json(fig, validate=False)]

        # Generate predictions for the entire dataset
        assignments = pd.DataFrame(dict(
            actual=pred,
            predicted=clf.predict(df),
            predict_log_proba=clf.predict_log_proba(df).tolist(),
            predict_proba=clf.predict_proba(df).tolist(),
            group=pd.Series(
                [
                    (
                        "training"
                        if ix in train_ix
                        else (
                            "testing"
                            if ix in test_ix
                            else None
                        )
                    )
                    for ix in df.index.values
                ],
                index=df.index
            )
        ))
        self.save_results("assignments", assignments, figures=figures)

        # Make a Series with the weights for each variable
        # in the model
        weights = pd.DataFrame(
            dict(
                weights=pd.Series(
                    clf.feature_importances_,
                    index=df.columns
                )
            )
        )
        self.save_results("weights", weights, figures=figures)
