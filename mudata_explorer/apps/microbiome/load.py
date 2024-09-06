from biom import load_table
import pandas as pd
from mudata_explorer.helpers import cirro
from mudata_explorer.parsers import curatedMetagenomicData
from mudata_explorer.parsers import microbiome
from mudata_explorer.apps.helpers import load_mdata
import streamlit as st
from tempfile import NamedTemporaryFile


def _load_data_biom():
    st.write('Load Data From BIOM File')

    biom = st.file_uploader("Upload a BIOM file", type=["biom"])
    if biom is None:
        return
    
    with NamedTemporaryFile() as f:
        f.write(biom.getvalue())
        f.seek(0)
        table = load_table(f.name)
    adata = table.to_anndata()

    mdata = microbiome.parse_adata(adata, groupby_var=True)
    load_mdata(mdata)


def _load_data_cmd():
    st.markdown(
        '''Load Curated Metagenomic Data
        
Analyze abundance data produced by the
[Curated Metagenomic Data](https://waldronlab.io/curatedMetagenomicData/)
project.''')

    file = st.file_uploader("Curated Metagenomic Data")
    if file is None:
        return

    # Parse the curatedMetagenomicData format as AnnData
    adata = curatedMetagenomicData.parse_df(pd.read_csv(file, sep="\t"))

    # Run the microbiome analysis
    mdata = microbiome.parse_adata(adata, groupby_var=True)
    load_mdata(mdata)


def run():
    source = st.selectbox(
        'Load Data From',
        options=[
            'BIOM File',
            'Curated Metagenomic Data',
            'Cirro'
        ],
        index=None
    )
    st.session_state['microbiome-empty-input'] = False
    if source == 'BIOM File':
        _load_data_biom()
    elif source == 'Curated Metagenomic Data':
        _load_data_cmd()
    elif source == 'Cirro':
        cirro.load_from_cirro(
            filter_process_ids=[
                "ingest_biom",
                "process-nf-core-ampliseq-2-4-0",
                "curated_metagenomic_data"
            ],
            show_link=False
        )
    else:
        st.session_state['microbiome-empty-input'] = True
