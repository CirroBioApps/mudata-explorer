from sklearn.decomposition import PCA
import pandas as pd
from mudata_explorer.base.process import Process


class RunPCA(Process):

    type = "pca"
    name = "PCA"
    desc = "Principle Coordinates Analysis (PCA)"
    category = "Dimensionality Reduction"
    schema = {
        "table": {
            "type": "object",
            "label": "Data Table",
            "properties": {
                "data": {
                    "type": "dataframe",
                    "select_columns": True,
                    "query": "",
                }
            }
        },
        "outputs": {
            "type": "object",
            "label": "Outputs",
            "properties": {
                "dest_key": {
                    "type": "string",
                    "default": "X_pca",
                    "label": "Label to use for results",
                    "help": """
                    Key to use when saving the output to the container
                    """
                }
            }
        }
    }
    outputs = {
        "coords": {
            "type": pd.DataFrame,
            "label": "PCA Coordinates",
            "desc": "Table of scalar values assigned to each component",
            "modality": "table.data.tables",
            "axis": "table.data.axis",
            "attr": "outputs.dest_key"
        }
    }

    def execute(self) -> pd.DataFrame:

        df: pd.DataFrame = self.params["table.data.dataframe"]

        # Run PCA
        pca = PCA()
        ndim = min(df.shape[0], df.shape[1])
        res = pd.DataFrame(
            pca.fit_transform(df),
            index=df.index,
            columns=[f"PC{i + 1}" for i in range(ndim)]
        )

        self.save_results("coords", res)
