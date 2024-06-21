import pandas as pd
from scipy import stats
from mudata_explorer.base.process import Process


class RunAnova(Process):

    type = "anova"
    name = "One-Way ANOVA"
    help_text = """
    The one-way ANOVA is used to determine if there are any statistically
    significant differences between the means of two or more independent
    (unrelated) groups. It assumes that the data is normally distributed
    and that the variances of the groups are equal.

    If these assumptions are not met, the results of the ANOVA may not be
    valid. In this case, consider using a non-parametric alternative such
    as the Kruskal-Wallis test.

    Wikipedia: [One-Way ANOVA](https://en.wikipedia.org/wiki/Analysis_of_variance)
    Scipy: [scipy.stats.f_oneway](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.f_oneway.html)

    Results will include:
    
    - F-statistic: The F-statistic of the ANOVA
    - pvalue: The p-value of the ANOVA
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
                    "help": "Select the data which will be used to perform the ANOVA",
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
                    "default": "anova",
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
            "label": "ANOVA Results",
            "desc": "Results of the one-way ANOVA test.",
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
                **self.run_anova(cvals, group)
            )
            for cname, cvals in df.items()
        ]).set_index("index")

        self.save_results("results", res)

    @staticmethod
    def run_anova(vals, group):

        res = stats.f_oneway(
            *[vals[group == g] for g in group.unique()]
        )
        return dict(
            f_statistic=res.statistic,
            pvalue=res.pvalue
        )
