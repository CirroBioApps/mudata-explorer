import pandas as pd
from scipy import stats
from mudata_explorer.base.process import Process


class RunKruskal(Process):

    type = "kruskal"
    name = "Kruskal-Wallis H-test"
    help_text = """
The Kruskal-Wallis H-test tests the null hypothesis that the population
median of all of the groups are equal. It is a non-parametric alternative
to the one-way ANOVA and extends the Mann-Whitney U test to more than two groups.

Documentation:

- Wikipedia: [Kruskal-Wallis test](https://en.wikipedia.org/wiki/Kruskal%E2%80%93Wallis_one-way_analysis_of_variance)
- Scipy: [scipy.stats.kruskal](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.kruskal.html)
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
                "grouping": {
                    "type": "dataframe",
                    "label": "Grouping",
                    "help": "Select the column containing categories to compare", # noqa
                    "columns": {"grouping": {"label": "Grouping"}}
                }
            }
        },
        "outputs": {
            "type": "object",
            "label": "Outputs",
            "properties": {
                "dest_key": {
                    "type": "string",
                    "default": "kruskal",
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
            "label": "Kruskal-Wallis Results",
            "desc": "Results of the one-way Kruskal-Wallis H-test.",
            "modality": "table.data.tables",
            "axis": "table.data.axis.T",
            "attr": "outputs.dest_key"
        }
    }

    def execute(self):

        df: pd.DataFrame = self.params["table.data.dataframe"]
        df = df.dropna()
        assert df.shape[0] > 0, "Null values in all rows."

        group: pd.Series = (
            self.params
            ["table.grouping.dataframe"]
            ["grouping"]
            .dropna()
        )

        # Get the set of observations which have both data and grouping
        shared = df.index.intersection(group.index)

        assert len(shared) > 0, "No shared rows between data and grouping"

        df = df.loc[shared]
        group = group.loc[shared]

        # Iterate over each variable in the data to run the ANOVA
        res = pd.DataFrame([
            dict(
                index=cname,
                **self.run_kruskal(cvals, group)
            )
            for cname, cvals in df.items()
            if cvals.nunique() > 1
        ]).set_index("index")

        # Assign a rank order to the results such that the
        # largest f-statistic is ranked first
        res["rank"] = res["f_statistic"].rank(ascending=False)

        self.save_results("results", res)

    @staticmethod
    def run_kruskal(vals, group):

        res = stats.kruskal(
            *[vals[group == g] for g in group.unique()]
        )
        return dict(
            f_statistic=res.statistic,
            pvalue=res.pvalue
        )
