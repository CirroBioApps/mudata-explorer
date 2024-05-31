from typing import Union
import pandas as pd
from sklearn.cluster import DBSCAN
from mudata_explorer.base.process import Process


class RunDBSCAN(Process):

    type = "dbscan"
    name = "DBSCAN"
    desc = "DBSCAN Clustering"
    categories = ["Clustering"]
    output_type = pd.Series
    schema = {
        "data": {
            "type": "dataframe",
            "select_columns": True,
            "query": "",
        },
        "use_zscore": {
            "type": "boolean",
            "label": "Use Z-Score Normalization",
            "help": """
            Normalize the data before clustering. This is useful when the
            features have different scales.
            """
        },
        "eps": {
            "type": "float",
            "min_value": 0.,
            "default": 0.5,
            "label": "Distance threshold - eps",
            "help": """
            The maximum distance between two samples for one to be considered
            as in the neighborhood of the other. This is not a maximum bound
            on the distances of points within a cluster. This is the most
            important DBSCAN parameter to choose appropriately for your data set
            and distance function.
            """
        },
        "min_samples": {
            "type": "integer",
            "min_value": 2,
            "default": 5,
            "label": "Neighborhood Size - min_samples",
            "help": """
            The number of samples (or total weight) in a neighborhood for a point to
            be considered as a core point. This includes the point itself. If
            `min_samples` is set to a higher value, DBSCAN will find denser clusters,
            whereas if it is set to a lower value, the found clusters will be more
            sparse.
            """
        },
        "metric": {
            "type": "string",
            "label": "DBSCAN: Metric",
            "enum": ["cosine", "euclidean", "manhattan", "correlation", "jaccard"],
            "help": """
            The metric to use when calculating distance between instances in a
            feature array.
            """
        }
    }

    def execute(self) -> Union[pd.Series, pd.DataFrame]:

        df: pd.DataFrame = self.params["data.dataframe"].dropna()
        msg = "Null values in all rows - remove invalid columns"
        assert df.shape[0] > 0, msg

        clusters = (
            DBSCAN(
                eps=self.params["eps"],
                min_samples=int(self.params["min_samples"]),
                metric=self.params["metric"]
            )
            .fit_predict(
                df.values
            )
        )
        return pd.Series(
            list(map(str, clusters)),
            index=df.index
        )
