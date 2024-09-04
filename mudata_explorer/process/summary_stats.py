import pandas as pd
from mudata_explorer.base.process import Process


class SummaryStats(Process):

    type = "summary-stats"
    name = "Summary Statistics"
    help_text = """
    Calculate a variety of summary statistics for the selected data.

    Note that only numerical data may be summarized in this way.

    - count: Number of non-null values
    - prop_valid: Proportion of non-null values
    - prop_valid_rank: Rank ordering of prop_valid (highest first)
    - nunique: Number of unique values
    - nunique_rank: Rank ordering of number of unique values (highest first)
    - median: Median value
    - median_rank: Rank ordering of median value (highest first)
    - mean: Mean value
    - mean_rank: Rank ordering of mean value (highest first)
    - std: Standard deviation
    - min: Minimum value
    - max: Maximum value
    - 25%: 25th percentile
    - 50%: 50th percentile
    - 75%: 75th percentile
    - prop_positive: Proportion of samples with positive values
    - prop_negative: Proportion of samples with negative values

    """
    category = "Summary Statistics"
    schema = {
        "table": {
            "type": "object",
            "properties": {
                "data": {
                    "label": "Data Table",
                    "type": "dataframe",
                    "select_columns": True,
                    "query": "",
                    "dropna": False
                }
            }
        },
        "outputs": {
            "type": "object",
            "label": "Outputs",
            "properties": {
                "dest_key": {
                    "type": "string",
                    "default": "summary_stats",
                    "label": "Label to use for results",
                    "help": """
                    Key to use when saving the output
                    """
                }
            }
        }
    }
    outputs = {
        "res": {
            "type": pd.DataFrame,
            "label": "Summary Statistics",
            "desc": "Summary statistics for each column",
            "modality": "table.data.tables",
            "axis": "table.data.axis",
            "attr": "outputs.dest_key"
        }
    }

    def execute(self):

        df: pd.DataFrame = self.params["table.data.dataframe"]

        assert df is not None, self.params
        # Calculate summary statistics
        res = df.apply(self.summary_stats, axis=1)

        # Add the rank-order of various values
        for kw in ['nunique', 'median', 'mean', 'prop_valid']:
            res[f'{kw}_rank'] = res[kw].rank(ascending=False)

        # Save the results and the figure
        self.save_results(
            "res",
            res
        )

    @staticmethod
    def summary_stats(vals: pd.Series):

        output = vals.dropna().describe()
        output['median'] = vals.dropna().median()
        output['prop_valid'] = (
            vals.dropna().shape[0] / vals.shape[0]
        )
        output['nunique'] = vals.dropna().nunique()
        output['prop_positive'] = (vals > 0).mean()
        output['prop_negative'] = (vals < 0).mean()
        return output
