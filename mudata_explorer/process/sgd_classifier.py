import numpy as np
import pandas as pd
import plotly.express as px
from plotly import io
from sklearn import linear_model
from mudata_explorer.process.ml_classifier import MLClassifier


class SGDClassifier(MLClassifier):

    type = "sgd_classifier"
    name = "SGD Classifier"
    help_text = """
Linear classifiers (SVM, logistic regression, etc.) with SGD training.

This estimator implements regularized linear models with stochastic
gradient descent (SGD) learning: the gradient of the loss is estimated
each sample at a time and the model is updated along the way with a
decreasing strength schedule (aka learning rate). SGD allows minibatch
(online/out-of-core) learning via the partial_fit method.
For best results using the default learning rate schedule, the data should
have zero mean and unit variance (i.e. z-score).

This analysis is used to predict _categorical_ variables.
    
The SGD Classifier is implemented using the
scikit-learn library, and more complete information can
be found on the [scikit-learn documentation](https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.SGDClassifier.html).

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
                    "help": "Select the column containing categorical values to predict", # noqa
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
                "loss": {
                    "type": "string",
                    "default": "hinge",
                    "label": "Loss Function",
                    "help": "The loss function to be used.",
                    "enum": [
                        "hinge",
                        "log_loss",
                        "modified_huber",
                        "squared_hinge",
                        "perceptron",
                        "huber",
                        "epsilon_insensitive",
                        "squared_epsilon_insensitive"
                    ]
                },
                "penalty": {
                    "type": "string",
                    "default": "l2",
                    "label": "Penalty",
                    "help": "The penalty (aka regularization term) to be used.",
                    "enum": ["l2", "l1", "elasticnet"]
                },
                "alpha": {
                    "type": "float",
                    "label": "Alpha: Regularization Strength",
                    "help": "Regularization improves the conditioning of the problem and reduces the variance of the estimates. Larger values specify stronger regularization.", # noqa
                    "default": 1.0,
                    "min_value": 0.0
                },
                "l1_ratio": {
                    "type": "float",
                    "label": "L1 Ratio",
                    "help": "The Elastic Net mixing parameter, with 0 <= l1_ratio <= 1. l1_ratio=0 corresponds to L2 penalty, l1_ratio=1 to L1. Only used if penalty='elasticnet'.", # noqa
                    "default": 0.15,
                    "min_value": 0.0,
                    "max_value": 1.0
                },
                "balance_groups": {
                    "type": "string",
                    "label": "Balance Groups",
                    "help": "Whether to select a random subset with equal numbers per class.", # noqa
                    "enum": ["Use All Data", "Balance Groups"],
                    "default": "Use All Data"
                },
                "fit_intercept": {
                    "type": "boolean",
                    "default": True,
                    "label": "Fit Intercept",
                    "help": "Whether to calculate the intercept for this model." # noqa
                },
                "tol": {
                    "type": "string",
                    "default": "0.001",
                    "label": "Tolerance",
                    "help": "Stopping criterion."
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
                    "default": "sgd_classifier",
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

        # Get the data and predictor
        # This will also balance the groups if necessary
        df, pred = self.get_data_and_predictor()

        train_prop = self.params["model_params.train_prop"]
        random_state = self.params["model_params.random_state"]
        params = dict(
            loss=self.params["model_params.loss"],
            penalty=self.params["model_params.penalty"],
            l1_ratio=self.params["model_params.l1_ratio"],
            alpha=self.params["model_params.alpha"],
            fit_intercept=self.params["model_params.fit_intercept"],
            tol=float(self.params["model_params.tol"]),
            random_state=random_state
        )

        # Get the training and testing sets
        train_ix = df.sample(frac=train_prop, random_state=random_state).index
        test_ix = df.index.difference(train_ix)

        assert train_ix.shape[0] > 0, "No training data - increase train_prop"
        assert test_ix.shape[0] > 0, "No testing data - decrease train_prop"

        # Run the SGDClassifier
        clf = linear_model.SGDClassifier(**params)
        clf = clf.fit(df.loc[train_ix], pred.loc[train_ix])

        # Compute the performance of the model by comparing the
        # accuracy of classification on the held out testing set
        # against the accuracy observed under random permutation
        test_predictions = clf.predict(df.loc[test_ix])
        score, null_dist, perc = self.score_by_permutation(
            test_predictions,
            pred.loc[test_ix].tolist()
        )
        print(f"Model Score: {score:.2f} (Above {int(100 * perc)}% of null distribution)") # noqa

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

        assignments = pd.DataFrame(dict(
            actual=pred,
            predicted=predicted,
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
            {
                (
                    "weight"
                    if len(clf.classes_) == 2
                    else
                    f"weight: {cls_name}"
                ): pd.Series(
                    clf.coef_[cls_ix],
                    index=df.columns
                )
                for cls_ix, cls_name in enumerate(clf.classes_[1:])
            }
        )

        # Make a barplot showing the weights of each feature
        figures = [
            io.to_json(
                px.bar(
                    cvals.loc[cvals.fillna(0) != 0].sort_values(),
                    orientation="h",
                    title=(
                        f"Feature Weights (n={(cvals.fillna(0) != 0).sum():,})"
                        if cname == "weight"
                        else f"Feature Weights: {cname.split(': ')[1]} (n={(cvals != 0).sum():,})"
                    ),
                    labels={"value": "Weight", "index": "Feature"}
                ),
                validate=False
            )
            for cname, cvals in weights.items()
        ]

        # Also save the absolute value of the weights
        weights = pd.concat([
            weights,
            weights.abs().rename(columns=lambda x: f"abs_{x}")
        ], axis=1)

        # Save the maximum absolute weight, as well as the rank order
        weights = weights.assign(
            max_abs_weight=weights.abs().max(axis=1),
            rank=weights.abs().mean(axis=1).rank(ascending=False)
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
