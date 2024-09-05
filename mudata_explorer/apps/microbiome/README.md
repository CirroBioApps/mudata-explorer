[![example dataset](https://github.com/CirroBioApps/mudata-explorer/blob/main/demo_data/curatedMetagenomicData/data/HMP_2019_ibdmdb/hmp-2019-ibdmdb-c80d1c5c40d77a8e.png?raw=true)](https://mudata-explorer.streamlit.app/?file=https://github.com/CirroBioApps/mudata-explorer/blob/main/demo_data/curatedMetagenomicData/data/HMP_2019_ibdmdb/hmp-2019-ibdmdb-c80d1c5c40d77a8e.h5mu)

[Open example dataset - HMP IBDMDB (2019)](https://mudata-explorer.streamlit.app/views?file=https://github.com/CirroBioApps/mudata-explorer/raw/main/demo_data/curatedMetagenomicData/data/HMP_2019_ibdmdb/hmp-2019-ibdmdb-c80d1c5c40d77a8e.h5mu)

### Input Data:

- Microbial Abundance Table (Samples by Microbes)
- Metadata Table (Information about Samples)

### Analysis Options:

- Taxonomic Level (Species, Genus, etc.)
- Filter Samples by Metadata
- Compare Communities by Metadata Annotation (Categorical or Continuous)

### Data Processing:

- UMAP Ordination
- Leiden Clustering
- Microbes Varying Across Clusters (Kruskal)
- Microbes Varying by Metadata (Kruskal for Categorical, Spearman for Continuous)
