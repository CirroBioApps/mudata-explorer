from mudata_explorer.base.process import Process
import pandas as pd


class MLClassifier(Process):
    """Base class for all machine learning classifiers."""

    def get_data_and_predictor(
        self,
        df_key="table.data.dataframe",
        pred_key="table.predictor.dataframe",
        pred_cname="predictor",
        balance_groups_key="model_params.balance_groups"
    ):

        df: pd.DataFrame = self.params[df_key]
        df = df.dropna()
        assert df.shape[0] > 0, "Null values in all rows."

        pred: pd.Series = (
            self.params
            [pred_key]
            [pred_cname]
            .dropna()
        )

        # Get the set of observations which have both data and predictor
        shared = df.index.intersection(pred.index)

        assert len(shared) > 0, "No shared rows between data and predictor"

        # If the user has selected to balance the groups
        if self.params[balance_groups_key] == "Balance Groups":
            pred = pred.loc[shared]
            # Get the counts of each class
            class_counts = pred.astype(str).value_counts()
            # Get the minimum number of samples in a class
            min_count = class_counts.min()
            assert min_count > 0, "Not enough samples in each class"
            # Randomly sample the same number of samples from each class
            shared = [
                *[
                    ix
                    for _, cat_pred in pred.groupby(pred)
                    if cat_pred.shape[0] >= min_count
                    for ix in cat_pred.sample(min_count).index
                ]
            ]

        return df.loc[shared], pred.loc[shared]
