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
    - nunique: Number of unique values
    - median: Median value
    - mean: Mean value
    - std: Standard deviation
    - min: Minimum value
    - max: Maximum value
    - 25%: 25th percentile
    - 50%: 50th percentile
    - 75%: 75th percentile

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

        # Calculate summary statistics
        res = df.apply(self.summary_stats, axis=1)

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
