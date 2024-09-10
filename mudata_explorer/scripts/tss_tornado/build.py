#!/usr/bin/env python

import json
from typing import List
import click
from pathlib import Path
import anndata as ad
from mudata import MuData
import pandas as pd
import numpy as np
import pyBigWig
from sklearn.cluster import KMeans
import umap
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("TSS Tornado")


def _read_data(
    file_paths: List[str],
    tss_coords: pd.DataFrame,
    suffix: str,
    window_size: int,
    nbins: int = 100
) -> ad.AnnData:

    # Make sure that each file ends with the suffix
    for fp in file_paths:
        assert fp.endswith(suffix), fp

    # Read in the coverage information for each .bigWig file
    # and return a merged AnnData object
    adata: ad.AnnData = ad.concat(
        [
            _calc_coverage(fp, tss_coords, window_size)
            for fp in file_paths
        ],
        label="sample",
        keys=[
            fp.split("/")[-1][:-len(suffix)]
            for fp in file_paths
        ]
    )
    adata.obs_names_make_unique()

    # Calculate the average depth across binned positions
    adata.obsm["binned_coverage"] = _calc_bins(adata, nbins)

    # Sort by 1D UMAP coordinates
    adata = adata[
        (
            _run_umap(
                adata.obsm["binned_coverage"],
                n_components=1
            )
            .sort_values("UMAP1")
            .index
        )
    ]

    # Add UMAP coordinates for the binned coverage
    adata.obsm["umap"] = _run_umap(adata.obsm["binned_coverage"])

    # Calculate the average depth at each position for each sample
    adata.uns['avg_coverage'] = (
        adata
        .to_df()
        .groupby(adata.obs["sample"])
        .mean()
        .T
    )

    # Run k-means clustering on the genes, across all samples
    adata.obsm['kmeans'] = pd.DataFrame({
        k: _run_kmeans(adata.to_df(), k)
        for k in range(2, 11)
    }, index=adata.obs_names).rename(columns=str)

    # Calculate the average depth per cluster, per sample, for each
    # value of K
    for k in range(2, 11):
        logger.info(f"Calculating Average Coverage per Cluster for K={k}")
        for cluster in adata.obsm['kmeans'][str(k)].unique():
            cluster_adata = adata[adata.obsm['kmeans'][str(k)] == cluster]
            adata.uns[f'avg_coverage.{k}.{cluster}'] = (
                cluster_adata.to_df()
                .groupby(cluster_adata.obs['sample'])
                .mean()
                .T
            )

    return adata


def _run_umap(df, n_components=2, **kwargs):
    return pd.DataFrame(
        umap.UMAP(n_components=n_components, **kwargs).fit_transform(df),
        columns=[f"UMAP{i + 1}" for i in range(n_components)],
        index=df.index
    )


def _calc_bins(adata: ad.AnnData, nbins: int):

    # Calculate which position goes into which bin
    pos_bins = pd.qcut(
        adata.var_names.astype(int),
        nbins,
        labels=False
    )

    # Take the average within each bin
    bin_avgs = (
        adata
        .to_df()
        .T
        .groupby(pos_bins)
        .mean()
        .T
    )
    assert bin_avgs.shape[1] == nbins, bin_avgs

    return bin_avgs


def _calc_coverage(
    fp: Path,
    tss_coords: pd.DataFrame,
    window_size: int,
    n_genes=100
) -> pd.DataFrame:

    # Read in the bigWig file
    logger.info(f"Reading BigWig File: {fp}")
    bw = pyBigWig.open(str(fp))

    # Subset the genes to those which have a TSS which is more than
    # window_size / 2 from the start or end of the chromosome

    half_window = int(window_size / 2)

    filtered_tss = tss_coords.loc[
        tss_coords.apply(
            lambda r: (
                bw.chroms(r['chrom']) is not None
                and (r["tss"] - half_window) >= 0
                and (r["tss"] + half_window) < bw.chroms(r["chrom"])
            ),
            axis=1
        )
    ]

    # Make a numeric matrix with coverage information
    # spanning the TSS for each gene
    logger.info(f"Getting coverage for {filtered_tss.shape[0]:,} genes")
    cov = pd.DataFrame(
        [
            bw.values(
                row["chrom"],
                row["tss"] - half_window,
                row["tss"] + half_window,
                numpy=True
            )
            for _, row in filtered_tss.iterrows()
        ],
        columns=range(-half_window, half_window),
        index=filtered_tss.index
    )

    # Fill in any NaN values with zeros
    cov = cov.fillna(0)

    # Take the top n_genes genes
    cov = cov.loc[
        (
            cov
            .sum(axis=1)
            .sort_values(ascending=False)
            .head(n_genes)
            .index
        )
    ]

    logger.info(f"Genes with Data: {cov.shape[0]:,} genes")

    # # Normalize to total sequencing depth
    # cov = window_size * cov / cov.sum().sum()

    # Normalize to sequencing depth per-gene
    cov = cov.apply(lambda x: x / x.sum(), axis=1)

    # Reformat as AnnData
    adata = ad.AnnData(X=cov, obs=filtered_tss.loc[cov.index])

    return adata


def _read_tss(bed_file: Path) -> pd.DataFrame:
    # Read the full BED file
    logger.info(f"Reading TSS Coordinates from {bed_file}")
    bed = pd.read_csv(
        bed_file,
        sep="\t",
        names=[
            "chrom", "start", "end", "name", "score", "strand",
            "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes",
            "blockStarts"
        ]
    )
    # Pick the TSS at the start or the end, depending on the strand
    bed = (
        bed
        .assign(
            tss=bed.apply(
                lambda r: r["start"] if r["strand"] == "+" else r["end"],
                axis=1
            )
        )
        .reindex(columns=["name", "chrom", "tss"])
    )
    logger.info(f"Read in TSS Coordinates: {bed.shape[0]:,} genes")

    # If there are any duplicate gene names
    if bed["name"].value_counts().max() > 1:
        bed = bed.groupby("name").first().reset_index()
        logger.info(f"Filtered to unique gene names: {bed.shape[0]:,} genes")

    # logger.info("Subsetting to 1000 genes for TESTING") # FIXME
    # return bed.head(1000)

    return bed


def _run_kmeans(df: pd.DataFrame, k: int):

    logger.info(f"Running K-Means Clustering (k={k})")

    clusters = (
        KMeans(n_clusters=k)
        .fit_predict(df.values)
    )
    return [
        f"Cluster {int(x) + 1}"
        for x in clusters
    ]


def _write_mdata(adata: ad.AnnData, output_prefix: str):

    # Format as a MuData object
    mdata = MuData(
        dict(
            binned_coverage=ad.AnnData(
                adata.obsm["binned_coverage"],
                obsm=dict(
                    umap=adata.obsm["umap"],
                    kmeans=adata.obsm["kmeans"]
                )
            )
        ),
        obs=adata.obs,
        uns=adata.uns
    )
    # Add the view type
    mdata.uns["mudata-explorer-views"] = json.dumps([
        dict(type="tss-tornado")
    ])

    fp = f"{output_prefix}.h5mu"
    mdata.write(fp)
    logger.info(f"MuData object written to {fp}")


@click.command()
@click.argument('file_paths', nargs=-1, type=click.Path(exists=True))
@click.option('--bed', type=click.Path(exists=True), help='Path to BED file')
@click.option('--suffix', type=str, help='File suffix', default='.bigWig')
@click.option('--window-size', type=int, help='Window Size', default=6000)
@click.option('--output', type=str, help='Output Prefix', default='output')
def main(
    file_paths,
    bed,
    suffix=".bigWig",
    window_size=6000,
    output="output"
):

    # Read in the TSS coordinates from a .bed file
    tss_coords = _read_tss(bed)

    # Read in coverage data for each dataset as an AnnData object
    # with coverage information spanning `window_size` around each TSS.
    # The total coverage will be normalized to total sequencing depth
    # on a per-sample basis.
    # The `sample` column will be used to differentiate samples.
    adata = _read_data(file_paths, tss_coords, suffix, window_size)

    # Write output as MuData
    _write_mdata(adata, output)


if __name__ == "__main__":
    main()
