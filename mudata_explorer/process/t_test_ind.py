import numpy as np
import pandas as pd
from scipy import stats
from mudata_explorer.base.process import Process


class TTestIndependent(Process):

    type = "t-test-ind"
    name = "T-test (Independent Samples)"
    help_text = """
    The independent samples t-test is used to determine if there is a significant
    difference between the means of two independent groups. It assumes that the
    data is normally distributed and that the variances of the two groups are equal.

    > Make sure that the 'Group' column only contains two unique values.

    Results will include:

    - pvalue: The p-value of the t-test
    - neg_log10_pvalue: The negative log10 of the p-value
    - t_statistic: The t-statistic of the t-test
    - median_diff: The difference in medians between the two groups

    """ # noqa
    category = "Statistical Tests"
    schema = {
        "data": {
            "type": "dataframe",
            "label": "Data",
            "help": "Select the data used to perform the t-test",
            "select_columns": True,
            "query": "",
            "dropna": False
        },
        "group": {
            "label": "Group",
            "help": "Select the column containing the groups to compare between",
            "type": "dataframe",
            "columns": {"group": {"label": "Group"}}
        },
        "dest_key": {
            "type": "string",
            "default": "t_test_ind",
            "label": "Label to use for results",
            "help": """
            Key to use when saving the output of the t-test.
            """
        }
    }
    outputs = {
        "res": {
            "type": pd.DataFrame,
            "label": "Analysis Results",
            "desc": "Results of the t-test.",
            "modality": "data.tables",
            "axis": "data.axis.T",
            "attr": "dest_key"
        }
    }

    def execute(self):

        df: pd.DataFrame = self.params["data.dataframe"]

        group: pd.Series = (
            self.params
            ["group.dataframe"]
            ["group"]
            .dropna()
        )

        # Get the set of observations which have both data and predictor
        shared = df.index.intersection(group.index)

        assert len(shared) > 0, "No shared rows between data and group"

        df = df.loc[shared]
        group = group.loc[shared]

        # The predictor must be a binary variable
        assert group.nunique() == 2, "Group must have two unique values"

        # Perform the t-test
        group1 = df[group == group.unique()[0]]
        group2 = df[group == group.unique()[1]]
        res = pd.DataFrame([
            ttest_ind(group1[col].dropna(), group2[col].dropna())
            for col in df.columns
        ], index=df.columns)

        # When calculating the -log10 p-value, protect
        # against zero values
        min_p = res.query("pvalue > 0")["pvalue"].min()
        res = res.assign(
            neg_log10_pvalue=(
                -1 * res["pvalue"].clip(lower=min_p).apply(np.log10)
            )
        )

        self.save_results("res", res)


def ttest_ind(a, b):
    res = stats.ttest_ind(a, b)
    return {
        "pvalue": res.pvalue,
        "t_statistic": res.statistic,
        "median_diff": a.median() - b.median()
    }
