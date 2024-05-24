from sklearn.cluster import DBSCAN
from mudata_explorer.base.process import Process


class RunDBSCAN(Process):

    type = "dbscan"
    name = "DBSCAN"
    desc = "DBSCAN Clustering"
    categories = ["Clustering"]
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
            "type": "number",
            "min_value": 0.,
            "value": 0.5,
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
            "type": "number",
            "min_value": 2,
            "step": 1,
            "value": 5,
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
        },
        "dest_key": {
            "type": "string",
            "label": "Destination Key",
            "help": "The name of the column which will be used for the results.",
            "value": "cluster"
        }
    }

    def execute(self):

        clusters = (
            DBSCAN(
                self.params["data.dataframe"],
                eps=self.params["eps"],
                min_samples=self.params["min_samples"],
                metric=self.params["metric"]
            )
            .fit_predict(
                self.params["data.dataframe"].values
            )
        )
        return list(map(str, clusters))
