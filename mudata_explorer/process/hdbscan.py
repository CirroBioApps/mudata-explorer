import pandas as pd
from sklearn.cluster import HDBSCAN
from mudata_explorer.base.process import Process


class RunHDBSCAN(Process):

    type = "hdbscan"
    name = "HDBSCAN Clustering"
    help_text = """
    Cluster data using hierarchical density-based clustering.

    [HDBSCAN - Hierarchical Density-Based Spatial Clustering of Applications with Noise](https://scikit-learn.org/stable/modules/generated/sklearn.cluster.HDBSCAN.html#sklearn.cluster.HDBSCAN).
    Performs DBSCAN over varying epsilon values and integrates the result to find a
    clustering that gives the best stability over epsilon. This allows HDBSCAN to
    find clusters of varying densities (unlike DBSCAN), and be more robust to
    parameter selection. Read more in the [User Guide](https://scikit-learn.org/stable/modules/clustering.html#hdbscan).

    For an example of how to use HDBSCAN, as well as a comparison to DBSCAN, please see the
    [plotting demo](https://scikit-learn.org/stable/auto_examples/cluster/plot_hdbscan.html#sphx-glr-auto-examples-cluster-plot-hdbscan-py).
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
                "min_cluster_size": {
                    "type": "integer",
                    "min_value": 1,
                    "default": 5,
                    "label": "Minimum Cluster Size",
                    "help": """
                    The minimum number of samples in a group for
                    that group to be considered a cluster; groupings
                    smaller than this size will be left as noise.
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
                "cluster_selection_epsilon": {
                    "type": "float",
                    "min_value": 0.,
                    "default": 0.,
                    "label": "Cluster Selection Epsilon",
                    "help": """
                    A distance threshold. Clusters below this value will be merged. See
                    [ref](https://scikit-learn.org/stable/modules/generated/sklearn.cluster.HDBSCAN.html#r6f313792b2b7-5)
                    for more information.
                    """ # noqa
                },
                "metric": {
                    "type": "string",
                    "label": "DBSCAN: Metric",
                    "default": "euclidean",
                    "enum": [
                        'braycurtis',
                        'canberra',
                        'chebyshev',
                        'cityblock'
                        'dice',
                        'euclidean',
                        'hamming',
                        'haversine',
                        'infinity',
                        'jaccard',
                        'l1',
                        'l2',
                        'manhattan',
                        'minkowski',
                        'p',
                        'pyfunc',
                        'rogerstanimoto',
                        'russellrao',
                        'seuclidean',
                        'sokalmichener',
                        'sokalsneath'
                    ],
                    "help": """
                    The metric to use when calculating distance between
                    instances in a feature array.
                    """
                },
                "alpha": {
                    "type": "float",
                    "min_value": 0.,
                    "default": 1.0,
                    "label": "Alpha",
                    "help": """
                    A distance scaling parameter. See
                    [ref](https://scikit-learn.org/stable/modules/generated/sklearn.cluster.HDBSCAN.html#r6f313792b2b7-3)
                    for more information.
                    """
                },
                "algorithm": {
                    "type": "string",
                    "label": "Algorithm",
                    "enum": [
                        "auto",
                        "brute",
                        "kd_tree",
                        "ball_tree"
                    ],
                    "default": "auto”",
                    "help": """
                    Exactly which algorithm to use for computing core
                    distances; By default this is set to "auto" which attempts
                    to use a `KDTree` tree if possible, otherwise it uses a
                    `BallTree` tree. Both "kd_tree" and "ball_tree" algorithms
                    use the NearestNeighbors estimator.
                    """
                },
                "leaf_size": {
                    "type": "integer",
                    "min_value": 1,
                    "default": 40,
                    "label": "Leaf Size",
                    "help": """
                    Leaf size for trees responsible for fast nearest neighbour
                    queries when a `KDTree` or a `BallTree` are used as
                    core-distance algorithms. A large dataset size and small
                    leaf_size may induce excessive memory usage. If you are
                    running out of memory consider increasing the `leaf_size`
                    parameter. Ignored for algorithm="brute".
                    """
                },
                "cluster_selection_method": {
                    "type": "string",
                    "label": "Cluster Selection Method",
                    "enum": [
                        "eom",
                        "leaf"
                    ],
                    "default": "eom",
                    "help": """
                    The method used to select clusters from the condensed tree.
                    """
                },
                "allow_single_cluster": {
                    "type": "boolean",
                    "default": False,
                    "label": "Allow Single Cluster",
                    "help": """
                    By default HDBSCAN requires that at least two clusters be
                    formed. If this parameter is set to True then single
                    cluster results will be allowed.
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

        df: pd.DataFrame = self.params["table.data.dataframe"].dropna()
        msg = "Null values in all rows - remove invalid columns"
        assert df.shape[0] > 0, msg

        clusters = (
            HDBSCAN(
                min_cluster_size=self.params["clustering.min_cluster_size"],
                min_samples=int(self.params["clustering.min_samples"]),
                cluster_selection_epsilon=self.params["clustering.cluster_selection_epsilon"],
                metric=self.params["clustering.metric"],
                alpha=self.params["clustering.alpha"],
                algorithm=self.params["clustering.algorithm"],
                leaf_size=self.params["clustering.leaf_size"],
                cluster_selection_method=self.params["clustering.cluster_selection_method"],
                allow_single_cluster=self.params["clustering.allow_single_cluster"]
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

        self.save_results("res", res)
