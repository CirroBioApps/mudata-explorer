from sklearn.decomposition import PCA
import pandas as pd
import streamlit as st
from streamlit.delta_generator import DeltaGenerator
from mudata_explorer import app
from mudata_explorer.base.process import Process


class RunPCA(Process):

    type = "pca"
    name = "PCA"
    desc = "Principle Coordinates Analysis (PCA)"
    category = "Dimensionality Reduction"
    output_type = pd.DataFrame
    schema = {
        "data": {
            "type": "dataframe",
            "select_columns": True,
            "query": "",
        }
    }

    def execute(self) -> pd.DataFrame:

        df: pd.DataFrame = self.params["data.dataframe"]
        # Run PCA
        pca = PCA()
        ndim = min(df.shape[0], df.shape[1])
        return pd.DataFrame(
            pca.fit_transform(df),
            index=df.index,
            columns=[f"PC{i + 1}" for i in range(ndim)]
        )
