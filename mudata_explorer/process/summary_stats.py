import pandas as pd
from mudata_explorer.base.process import Process


class SummaryStats(Process):

    type = "summary-stats"
    name = "Summary Statistics"
    help_text = """
    Calculate a variety of summary statistics for the selected data.

    Note that only numerical columns may be summarized in this way.

    - count: Number of non-null values in each column
    - prop_valid: Proportion of non-null values in each column
    - nunique: Number of unique values in each column
    - median: Median value of each column
    - mean: Mean value of each column
    - std: Standard deviation of each column
    - min: Minimum value of each column
    - max: Maximum value of each column
    - 25%: 25th percentile of each column
    - 50%: 50th percentile of each column
    - 75%: 75th percentile of each column

    """ # noqa
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
            "axis": "table.data.axis.T",
            "attr": "outputs.dest_key"
        }
    }

    def execute(self):

        df: pd.DataFrame = self.params["table.data.dataframe"]

        # Calculate summary statistics
        res = df.apply(self.summary_stats, axis=0).T

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
        return output
