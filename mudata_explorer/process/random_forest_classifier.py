import numpy as np
import pandas as pd
import plotly.express as px
from plotly import io
from sklearn import ensemble
from mudata_explorer.base.process import Process


class RandomForestClassifier(Process):

    type = "random_forest_classifier"
    name = "Random Forest Classifier"
    help_text = """
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
                    "type": "dataframe",
                    "label": "Predictor",
                    "help": "Select the column containing categorical values to predict",
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
                "n_estimators": {
                    "type": "integer",
                    "label": "Number of Estimators",
                    "help": "The number of trees in the forest.",
                    "default": 100
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
                    "default": "random_forest_classifier",
                    "label": "Label to use for results",
                    "help": """
                    Key to use when saving the output to the container
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
        n_estimators = self.params["model_params.n_estimators"]
        criterion = self.params["model_params.criterion"]
        max_depth = self.params["model_params.max_depth"]
        min_samples_split = self.params["model_params.min_samples_split"]
        random_state = self.params["model_params.random_state"]

        # Get the training and testing sets
        train_ix = df.sample(frac=train_prop, random_state=random_state).index
        test_ix = df.index.difference(train_ix)

        assert train_ix.shape[0] > 0, "No training data - increase train_prop"
        assert test_ix.shape[0] > 0, "No testing data - decrease train_prop"

        # Run the RandomForestClassifier
        clf = ensemble.RandomForestClassifier(
            n_estimators=n_estimators,
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
        test_predictions = clf.predict(df.loc[test_ix])
        score, null_dist, perc = self.score_by_permutation(
            test_predictions,
            pred.loc[test_ix].tolist()
        )
        print(f"Model Score: {score:.2f} (Above {int(100 * perc)}% of null distribution)")

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
        fig.add_annotation(
            x=score,
            xanchor="left",
            y=-20,
            text=f"Model Score: {score:.2f} (Above {int(100 * perc)}% of null distribution)", # noqa
            showarrow=False
        )
        figures = [io.to_json(fig, validate=False)]

        # Generate predictions for the entire dataset
        predicted = clf.predict(df)
        # Get the probability of the predicted value
        predict_proba = [
            proba[list(clf.classes_).index(pred)]
            for pred, proba in zip(
                predicted,
                clf.predict_proba(df).tolist()
            )
        ]
        assignments = pd.DataFrame(dict(
            actual=pred,
            predicted=predicted,
            predict_proba=predict_proba,
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

    def score_by_permutation(self, pred, truth, n_reps=1000):
        null_dist = []
        for i in range(n_reps):
            null_dist.append(
                self.score(
                    pred,
                    truth,
                    permute=True
                )
            )
        score = self.score(pred, truth)
        perc = sum([n < score for n in null_dist]) / n_reps

        return score, null_dist, perc

    @staticmethod
    def score(pred, truth, permute=False):

        return np.mean([
            x == y
            for x, y in zip(
                pred,
                (
                    pd.Series(truth).sample(frac=1).tolist()
                    if permute
                    else truth
                )
            )
        ])
