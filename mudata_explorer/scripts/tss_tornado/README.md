# Transcription Start Site (TSS) Tornado Builder

Purpose: Build the data assets needed to run the TSS Tornado visualization.

### Data Sources

- Coverage data in .bigWig format (one per sample)
- Gene coordinates in .bed format

### Output Format

The output is a MuData object with the following structure:

```
  obs:  'name', 'chrom', 'tss', 'sample'
  uns:  'avg_coverage', 'avg_coverage.2.Cluster 1', 'avg_coverage.2.Cluster 2', 'avg_coverage.3.Cluster 2', 'avg_coverage.3.Cluster 3', 'avg_coverage.3.Cluster 1', 'avg_coverage.4.Cluster 4', 'avg_coverage.4.Cluster 3', 'avg_coverage.4.Cluster 1', 'avg_coverage.4.Cluster 2', 'avg_coverage.5.Cluster 4', 'avg_coverage.5.Cluster 5', 'avg_coverage.5.Cluster 1', 'avg_coverage.5.Cluster 3', 'avg_coverage.5.Cluster 2', 'avg_coverage.6.Cluster 1', 'avg_coverage.6.Cluster 6', 'avg_coverage.6.Cluster 3', 'avg_coverage.6.Cluster 5', 'avg_coverage.6.Cluster 2', 'avg_coverage.6.Cluster 4', 'avg_coverage.7.Cluster 1', 'avg_coverage.7.Cluster 5', 'avg_coverage.7.Cluster 6', 'avg_coverage.7.Cluster 4', 'avg_coverage.7.Cluster 2', 'avg_coverage.7.Cluster 7', 'avg_coverage.7.Cluster 3', 'avg_coverage.8.Cluster 2', 'avg_coverage.8.Cluster 5', 'avg_coverage.8.Cluster 7', 'avg_coverage.8.Cluster 8', 'avg_coverage.8.Cluster 1', 'avg_coverage.8.Cluster 6', 'avg_coverage.8.Cluster 3', 'avg_coverage.8.Cluster 4', 'avg_coverage.9.Cluster 7', 'avg_coverage.9.Cluster 9', 'avg_coverage.9.Cluster 8', 'avg_coverage.9.Cluster 5', 'avg_coverage.9.Cluster 4', 'avg_coverage.9.Cluster 1', 'avg_coverage.9.Cluster 3', 'avg_coverage.9.Cluster 6', 'avg_coverage.9.Cluster 2', 'avg_coverage.10.Cluster 8', 'avg_coverage.10.Cluster 9', 'avg_coverage.10.Cluster 10', 'avg_coverage.10.Cluster 2', 'avg_coverage.10.Cluster 7', 'avg_coverage.10.Cluster 1', 'avg_coverage.10.Cluster 5', 'avg_coverage.10.Cluster 6', 'avg_coverage.10.Cluster 4', 'avg_coverage.10.Cluster 3', 'mudata-explorer-views'
  1 modality
    binned_coverage:    9791 x 100
      obsm:     'umap', 'kmeans'
```