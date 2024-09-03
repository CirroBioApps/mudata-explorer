import pandas as pd
import scanpy as sc
from mudata_explorer.base.process import Process


class RunLeiden(Process):

    type = "leiden"
    name = "Leiden Clustering"
    help_text = """
    Leiden clustering is a graph-based clustering algorithm that is
    particularly well-suited to single-cell RNA-seq data. It is an
    improved version of the Louvain algorithm that is more robust to
    resolution parameters and can find clusters of varying sizes.

    The primary parameter used to tune the size of clusters generated
    by Leiden is the `resolution` parameter. Higher values of `resolution`
    will lead to more and smaller clusters, while lower values will lead
    to fewer and larger clusters.

    - [Traag, Waltman, and van Eck, 2019](https://www.nature.com/articles/s41598-019-41695-z)
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
                "resolution": {
                    "type": "float",
                    "min_value": 0.,
                    "default": 1.0,
                    "label": "Resolution",
                    "help": """
                    A parameter value controlling the coarseness of the
                    clustering. Higher values lead to more clusters.
                    """
                },
                "n_neighbors": {
                    "type": "integer",
                    "min_value": 2,
                    "default": 15,
                    "label": "Number of Neighbors",
                    "help": """
                    The size of local neighborhood (in terms of number of neighboring data points) used for manifold approximation. Larger values result in more global views of the manifold, while smaller values result in more local data being preserved. In general values should be in the range 2 to 100.
                    """ # noqa
                },
                "metric": {
                    "type": "string",
                    "label": "Distance Metric",
                    "default": "euclidean",
                    "enum": [
                        "braycurtis",
                        "canberra",
                        "chebyshev",
                        "cityblock",
                        "correlation",
                        "cosine",
                        "dice",
                        "euclidean",
                        "hamming",
                        "jaccard",
                        "kulsinski",
                        "l1",
                        "l2",
                        "mahalanobis",
                        "manhattan",
                        "minkowski",
                        "rogerstanimoto",
                        "russellrao",
                        "seuclidean",
                        "sokalmichener",
                        "sokalsneath",
                        "sqeuclidean",
                        "yule"
                    ],
                    "help": """
                    The metric to use when calculating distance between points.
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
                    "default": "leiden",
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

        # Build an AnnData object with the selected data
        adata = sc.AnnData(df)

        # Run the neighbors command using the selected metric
        sc.pp.neighbors(
            adata,
            n_neighbors=int(self.params["clustering.n_neighbors"]),
            metric=self.params["clustering.metric"],
            use_rep='X'
        )

        # Run Leiden clustering
        sc.tl.leiden(
            adata,
            resolution=self.params["clustering.resolution"],
            flavor="igraph",
            n_iterations=2,
            directed=False
        )

        # Format the results as a pandas series
        res = (
            adata
            .obs["leiden"]
            .apply(
                lambda x: f"Cluster {int(x) + 1}"
            )
        )

        self.save_results("res", res)
