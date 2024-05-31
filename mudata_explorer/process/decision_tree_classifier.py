import pandas as pd
import umap
from sklearn import tree
import streamlit as st
from mudata_explorer.base.process import Process


class DecisionTreeClassifier(Process):

    type = "decision_tree_classifier"
    name = "Decision Tree Classifier"
    desc = "Predict the value of a target variable by learning simple decision rules inferred from the data features"
    category = "Machine Learning"
    output_type = pd.DataFrame
    schema = {
        "data": {
            "label": "Data",
            "help": "Select the data used to train and test the model",
            "type": "dataframe",
            "select_columns": True,
            "query": "",
        },
        "predictor": {
            "label": "Predictor",
            "help": "Select the column to predict",
            "type": "dataframe",
            "columns": {"predictor": {"label": "Predictor"}}
        },
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
            "help": "The minimum number of samples required to split an internal node."
        },
        "random_state": {
            "type": "integer",
            "default": 0,
            "label": "Random State",
            "help": "Controls the randomness of the estimator."
        }
    }

    def execute(self) -> pd.DataFrame:

        df: pd.DataFrame = self.params["data.dataframe"]
        df = df.dropna()
        assert df.shape[0] > 0, "Null values in all rows."

        pred: pd.Series = (
            self.params
            ["predictor.dataframe"]
            ["predictor"]
            .dropna()
        )

        # Get the set of observations which have both data and predictor
        shared = df.index.intersection(pred.index)

        assert len(shared) > 0, "No shared rows between data and predictor"

        df = df.loc[shared]
        pred = pred.loc[shared]

        # Get the training and testing sets
        train = df.sample(
            frac=self.params["train_prop"],
            random_state=self.params["random_state"]
        )

        assert train.shape[0] > 0, "No training data - increase train_prop"

        # Run the DecisionTreeClassifier
        clf = tree.DecisionTreeClassifier(
            criterion=self.params["criterion"],
            max_depth=(
                self.params["max_depth"]
                if self.params["max_depth"] > 0
                else None
            ),
            min_samples_split=self.params["min_samples_split"],
            random_state=self.params["random_state"]
        )
        clf = clf.fit(train, pred.loc[train.index])

        # Generate predictions for the entire dataset
        output = pd.DataFrame(dict(
            actual=pred,
            predicted=clf.predict(df),
            group=pd.Series(
                [
                    "training" if ix in train.index else "testing"
                    for ix in df.index.values
                ],
                index=df.index
            )
        ))

        return output
