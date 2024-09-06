[![example dataset](https://github.com/CirroBioApps/mudata-explorer/blob/main/demo_data/curatedMetagenomicData/data/HMP_2019_ibdmdb/hmp-2019-ibdmdb-c80d1c5c40d77a8e.png?raw=true)](https://mudata-explorer.streamlit.app/views?file=https://github.com/CirroBioApps/mudata-explorer/raw/main/demo_data/curatedMetagenomicData/data/HMP_2019_ibdmdb/hmp-2019-ibdmdb-c80d1c5c40d77a8e.h5mu)

[Open example dataset - HMP IBDMDB (2019)](https://mudata-explorer.streamlit.app/views?file=https://github.com/CirroBioApps/mudata-explorer/raw/main/demo_data/curatedMetagenomicData/data/HMP_2019_ibdmdb/hmp-2019-ibdmdb-c80d1c5c40d77a8e.h5mu)

### Required Inputs:

- Microbial Abundance Table (Samples by Microbes)
- Metadata Table (Information about Samples)

### User Input:

- Select a taxonomic level for analysis (species, genus, etc.) (if available)
- Select which samples to include in the analysis
- Indicate groups of samples to compare, either by categorical (e.g. treatment vs. control) or continuous (e.g. age) variable

### Data Processing:

- Alpha Diversity (Shannon)
- Beta Diversity (Bray Curtis)
- Ordination (PCA, t-SNE, and UMAP)
- Unsupervised Clustering (Leiden)
- Identify microbes which vary across clusters
- Identify microbes which vary between groups, or by continuous variable
