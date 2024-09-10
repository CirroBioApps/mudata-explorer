from mudata_explorer.parsers import curatedMetagenomicData
from mudata_explorer.parsers import microbiome
from mudata_explorer.public_data.curatedMetagenomicData.data import datasets
from mudata_explorer.apps.helpers import load_mdata
from pathlib import Path
import pandas as pd
import streamlit as st


def parse(file: Path):
    
    # Parse the curatedMetagenomicData format as AnnData
    adata = curatedMetagenomicData.parse_df(pd.read_csv(file, sep="\t"))

    # Run the microbiome analysis
    mdata = microbiome.parse_adata(adata, groupby_var=True)
    load_mdata(mdata)


def run():
    
    # Let the user select a dataset
    dataset = st.selectbox(
        "Published Datasets",
        sorted(list(datasets.keys())),
        placeholder="Select a dataset",
        index=None
    )

    if dataset is not None:
        # Parse the selected dataset
        parse(datasets[dataset])
