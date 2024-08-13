import pandas as pd
from numpy import log10
from scipy import stats
from mudata_explorer.base.process import Process


class RunSpearman(Process):

    type = "spearman"
    name = "Spearman Rank Correlation"
    help_text = """
The Spearman rank-order correlation coefficient is a nonparametric measure of
the strength and direction of association between two ranked variables.
It assesses how well the relationship between two variables can be described
using a monotonic function.
A simple way of describing the analysis is that it is a Pearson correlation
on the ranks of the data (and so it does not assume a linear relationship).

Documentation:
- Wikipedia: [Spearman's rank correlation coefficient](https://en.wikipedia.org/wiki/Spearman%27s_rank_correlation_coefficient)
- Scipy: [scipy.stats.spearmanr](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.spearmanr.html)
    """ # noqa

    category = "Statistical Tests"
    schema = {
        "table": {
            "type": "object",
            "label": "Data Table",
            "properties": {
                "data": {
                    "type": "dataframe",
                    "label": "Data",
                    "help": "Select the data which will be used to perform the Kruskal-Wallis test", # noqa
                    "select_columns": True,
                    "query": "",
                },
                "comparitor": {
                    "type": "dataframe",
                    "label": "Comparitor",
                    "help": "Select the column containing values to compare", # noqa
                    "columns": {"comparitor": {"label": "Comparitor"}}
                }
            }
        },
        "outputs": {
            "type": "object",
            "label": "Outputs",
            "properties": {
                "dest_key": {
                    "type": "string",
                    "default": "spearman",
                    "label": "Label to use for results",
                    "help": """
                    Key to use when saving the output
                    """
                }
            }
        }
    }
    outputs = {
        "results": {
            "type": pd.DataFrame,
            "label": "Spearman Results",
            "desc": "Results of the Spearman rank sum correlation.",
            "modality": "table.data.tables",
            "axis": "table.data.axis.T",
            "attr": "outputs.dest_key"
        }
    }

    def execute(self):

        df: pd.DataFrame = self.params["table.data.dataframe"]
        df = df.dropna()
        assert df.shape[0] > 0, "Null values in all rows."

        comparitor: pd.Series = (
            self.params
            ["table.comparitor.dataframe"]
            ["comparitor"]
            .dropna()
        )

        # Get the set of observations which have both data and comparitor
        shared = df.index.intersection(comparitor.index)

        assert len(shared) > 0, "No shared rows between data and grouping"

        df = df.loc[shared]
        comparitor = comparitor.loc[shared]

        # Iterate over each variable in the data to run the test
        res = pd.DataFrame([
            dict(
                index=cname,
                mean=cvals.mean(),
                **self.run_spearman(cvals, comparitor)
            )
            for cname, cvals in df.items()
            if cvals.nunique() > 1
        ]).set_index("index")

        # Assign a rank order to the results such that the
        # largest f-statistic is ranked first
        res["statistic_rank"] = res["statistic"].rank(ascending=False)

        # Rank order by pvalue as well
        res["pvalue_rank"] = res["pvalue"].rank()

        # Calculate the -log10 pvalue
        res["neg_log10_pvalue"] = -res["pvalue"].apply(lambda x: log10(x))

        self.save_results("results", res)

    @staticmethod
    def run_spearman(vals, group):

        try:
            res = stats.spearmanr(
                vals,
                group
            )
        except ValueError as e:
            print(vals)
            raise e
        return dict(
            statistic=res.statistic,
            pvalue=res.pvalue
        )
