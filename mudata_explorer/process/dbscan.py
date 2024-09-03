import pandas as pd
from sklearn.cluster import DBSCAN
from mudata_explorer.base.process import Process


class RunDBSCAN(Process):

    type = "dbscan"
    name = "DBSCAN Clustering"
    help_text = """
    DBSCAN (Density-Based Spatial Clustering of Applications with Noise)
    is a clustering algorithm that groups together points that are
    closely packed together (points with many nearby neighbors), marking
    as outliers points that lie alone in low-density regions.

    The algorithm works by defining neighborhoods around each point and
    grouping points that are within a certain distance of each other.
    Points that are within the neighborhood of a core point are considered
    part of the same cluster. Points that are within the neighborhood of
    a non-core point but are not core points themselves are considered
    border points. Points that are not within the neighborhood of any
    core points are considered outliers.

    - [Wikipedia: DBSCAN](https://en.wikipedia.org/wiki/DBSCAN)
    - [Explanation of DBSCAN (video)](https://www.youtube.com/watch?v=C3r7tGRe2eI)
    """ # noqa
    category = "Clustering"
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
        "clustering": {
            "type": "object",
            "label": "Clustering Parameters",
            "properties": {
                "eps": {
                    "type": "float",
                    "min_value": 0.,
                    "default": 0.5,
                    "label": "Distance threshold - eps",
                    "help": """
                    The maximum distance between two samples for one to be
                    considered as in the neighborhood of the other.
                    This is not a maximum bound on the distances of points
                    within a cluster. This is the most important DBSCAN
                    parameter to choose appropriately for your data set
                    and distance function.
                    """
                },
                "min_samples": {
                    "type": "integer",
                    "min_value": 2,
                    "default": 5,
                    "label": "Neighborhood Size - min_samples",
                    "help": """
                    The number of samples (or total weight) in a neighborhood
                    for a point to be considered as a core point. This includes
                    the point itself. If `min_samples` is set to a higher
                    value, DBSCAN will find denser clusters, whereas if it is
                    set to a lower value, the found clusters will be more
                    sparse.
                    """
                },
                "metric": {
                    "type": "string",
                    "label": "DBSCAN: Metric",
                    "enum": [
                        "cosine",
                        "euclidean",
                        "manhattan",
                        "correlation",
                        "jaccard"
                    ],
                    "help": """
                    The metric to use when calculating distance between
                    instances in a feature array.
                    """
                }
            }
        },
        "outputs": {
            "type": "object",
            "label": "Outputs",
            "properties": {
                "dest_key": {
                    "type": "string",
                    "default": "dbscan",
                    "label": "Label to use for results",
                    "help": """
                    Key to use when saving the output of the clustering.
                    """
                }
            }
        }
    }
    outputs = {
        "res": {
            "type": pd.Series,
            "label": "Clusters",
            "desc": "Cluster assignments for each element",
            "modality": "table.data.tables",
            "axis": "table.data.axis",
            "attr": "outputs.dest_key"
        }
    }

    def execute(self):

        if self.params.get("table.data.dataframe") is None:
            raise Exception("Must select input data table")
        df: pd.DataFrame = self.params["table.data.dataframe"].dropna()
        msg = "Null values in all rows - remove invalid columns"
        assert df.shape[0] > 0, msg

        clusters = (
            DBSCAN(
                eps=self.params["clustering.eps"],
                min_samples=int(self.params["clustering.min_samples"]),
                metric=self.params["clustering.metric"]
            )
            .fit_predict(
                df.values
            )
        )
        res = pd.Series(
            list(map(str, clusters)),
            index=df.index
        )

        # Anything with a cluster assignment of -1 was unassigned
        res = res.loc[res != '-1']

        res = res.apply(
            lambda x: f"Cluster {int(x) + 1}"
        )

        self.save_results("res", res)
