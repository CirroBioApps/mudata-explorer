import pandas as pd
import numpy as np
from mudata_explorer.base.process import Process


class ShannonDiversity(Process):

    type = "shannon-diversity"
    name = "Shannon Diversity Index"
    help_text = """
Calculate the Shannon diversity index for the selected data.

The Shannon diversity index is a metric used in ecology to
estimate species diversity.
It is based on Claude Shannon's formula for entropy and takes
into account the number of species living in a habitat (richness)
and their relative abundance (evenness).
The index is denoted as H and is calculated as H = -Σpi * ln (pi).
It is typically used in environmental science to determine the species
biodiversity of a community.

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
                    "default": "shannon",
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
            "type": pd.Series,
            "label": "Shannon Diversity Index",
            "desc": "Shannon diversity index for each sample",
            "modality": "table.data.tables",
            "axis": "table.data.axis",
            "attr": "outputs.dest_key"
        }
    }

    def execute(self):

        df: pd.DataFrame = self.params["table.data.dataframe"]

        assert df is not None, self.params
        # Calculate summary statistics
        res = df.apply(self.shannon, axis=1)

        # Save the results and the figure
        self.save_results(
            "res",
            res
        )

    @staticmethod
    def shannon(vals: pd.Series):
        """Calculate the shannon index for a single sample."""
        vals = vals.dropna()
        total = vals.sum()
        if total == 0:
            return 0
        vals = vals[vals > 0]
        p: pd.Series = vals / total
        return -1 * (p * p.apply(np.log)).sum()
