from sklearn.decomposition import PCA
import pandas as pd
from mudata_explorer.base.process import Process


class RunPCA(Process):

    type = "pca"
    name = "PCA"
    desc = "Principle Coordinates Analysis (PCA)"
    help_text = """
    PCA is a linear dimensionality reduction technique that seeks to find the
    directions (or principal components) in which the data varies the most.
    These directions are found by computing the eigenvectors of the covariance
    matrix of the data.

    The first principal component is the direction in which the data varies
    the most. The second principal component is the direction in which the data
    varies the second most, and so on.

    The degree to which each principal component captures the variance in the
    data is reported in the results as a percentage.

    - [Wikipedia: Principal Component Analysis](https://en.wikipedia.org/wiki/Principal_component_analysis)
    - [Visual Explanation of PCA](https://towardsdatascience.com/principal-component-analysis-pca-explained-visually-with-zero-math-1cbf392b9e7d)
    """
    category = "Dimensionality Reduction"
    schema = {
        "table": {
            "type": "object",
            "properties": {
                "data": {
                    "label": "Data Table",
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

        # Make a new dataframe with the PCA results
        res = pd.DataFrame(
            pca.fit_transform(df),
            index=df.index,
            columns=[
                f"PC{i + 1} ({100 * pct:.2f}%)"
                for i, pct in enumerate(pca.explained_variance_ratio_)
            ]
        )

        self.save_results("coords", res)
