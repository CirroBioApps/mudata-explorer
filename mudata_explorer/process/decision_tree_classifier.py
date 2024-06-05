import pandas as pd
from sklearn import tree
from mudata_explorer.base.process import Process


class DecisionTreeClassifier(Process):

    type = "decision_tree_classifier"
    name = "Decision Tree Classifier"
    desc = """Predict the value of a target variable by learning
    simple decision rules inferred from the data features"""
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
            "type": pd.Series,
            "label": "Feature Weights",
            "desc": "Weights of each feature in the model.",
            "modality": "table.data.tables",
            "axis": "table.data.axis.T",
            "attr": "outputs.dest_key"
        }
    }

    def execute(self) -> pd.DataFrame:

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
        train = df.sample(frac=train_prop, random_state=random_state)

        assert train.shape[0] > 0, "No training data - increase train_prop"

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
        clf = clf.fit(train, pred.loc[train.index])

        # Generate predictions for the entire dataset
        assignments = pd.DataFrame(dict(
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
        self.save_results("assignments", assignments)

        # Make a Series with the weights for each variable
        # in the model
        weights = pd.Series(
            clf.feature_importances_,
            index=df.columns
        )
        self.save_results("weights", weights)
