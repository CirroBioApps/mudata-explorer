#!/usr/bin/env streamlit run

from anndata import AnnData
from sklearn.metrics import silhouette_score
import streamlit as st
from muon import MuData
import pandas as pd
import umap
import seaborn as sns
from sklearn import cluster
from typing import Tuple
import plotly.express as px

st.session_state["n_inputs"] = 5


def main():
    setup_inputs()
    setup_outputs()
    read_files()
    format_mudata()
    plot_clusters()


def plot_clusters():
    # Get the data to plot
    df, metadata = get_clusters_data()
    if df is None or metadata is None:
        tabs["cluster"].write("No data to plot")
        return

    # Plot a clustermap of the data
    g = sns.clustermap(df, cmap="Blues")
    tabs["cluster"].pyplot(g)

    # Cluster the data
    clusters = {
        n: run_clustering(df, n_clusters=n)
        for n in range(2, 11)
        if n < df.shape[0]
    }

    # Compute silhouette scores for the clusters
    silhouette_scores = {
        n: silhouette_score(df, clust)
        for n, clust in clusters.items()
    }

    # Plot the silhouette scores
    fig = px.line(
        x=list(silhouette_scores.keys()),
        y=list(silhouette_scores.values()),
        labels={"x": "Number of clusters", "y": "Silhouette score (1=best)"},
        title="Evaluate Clustering Performance"
    )
    tabs["cluster"].plotly_chart(fig)

    # Let the user select the number of clusters
    selected_n = tabs["cluster"].selectbox(
        "Select number of clusters",
        list(silhouette_scores.keys())
    )

    umap_coords = run_umap(
        df,
        n_neighbors=tabs["cluster"].number_input(
            "UMAP: Number of neighbors",
            value=15
        ),
        min_dist=tabs["cluster"].number_input(
            "UMAP: Minimum distance",
            value=0.1
        ),
        metric=tabs["cluster"].selectbox(
            "UMAP: Metric",
            ["cosine", "euclidean", "manhattan", "correlation", "jaccard"]
        )
    )

    # Merge the UMAP coordinates with the metadata
    plot_df = (
        umap_coords
        .assign(
            cluster=clusters[selected_n]
        )
        .merge(metadata, left_index=True, right_index=True)
    )

    # Plot the UMAP coordinates
    fig = px.scatter(
        data_frame=plot_df,
        x="UMAP1",
        y="UMAP2",
        color="cluster",
        hover_data=metadata.columns,
    )
    tabs["cluster"].plotly_chart(fig)

    # Summarize the selected markers for each cluster
    fig = px.box(
        (
            df
            .assign(cluster=clusters[selected_n])
            .melt(id_vars="cluster")
        ),
        color="cluster",
        y="value",
        x="variable",
    )
    tabs["cluster"].plotly_chart(fig)

    return

    # Let the user select a metadata column
    cnames = [
        cname for cname, cvals in metadata.items()
        if cvals.nunique() < 10 and cvals.nunique() > 1
    ]
    if len(cnames) > 0:
        selected_col = tabs["cluster"].selectbox(
            "Select metadata column to summarize",
            cnames
        )

        # Plot the proportion of observations in each cluster,
        # broken down by the selected metadata column
        df = plot_df.pivot_table(
            index=selected_col,
        )
        fig = px.histogram(
            plot_df,
            x=selected_col,
            color="cluster",
            barmode="overlay"
        )
        tabs["cluster"].plotly_chart(fig)


def get_clusters_data() -> Tuple[pd.DataFrame, pd.DataFrame]:
    mdata: MuData = st.session_state.get("mudata")
    if mdata is None:
        return None, None

    # Get the metadata table
    metadata = st.session_state.get("metadata_file_data")
    if metadata is None:
        return None, None

    # Let the user select the modality
    mod = tabs["cluster"].selectbox(
        "Select modality",
        list(mdata.mod.keys())
    )

    # Get the data for the selected modality
    df: pd.DataFrame = mdata.mod[mod].to_df()

    # Drop any rows which are entirely missing data
    df = df.reindex(index=[
        idx for idx, row in df.iterrows()
        if row.notnull().any()
    ])

    # Select the columns to include
    col_labels, default_col_labels = format_col_labels(df)

    # Let the user select the columns to include
    columns = tabs["cluster"].multiselect(
        "Select columns",
        col_labels,
        default=default_col_labels
    )
    # Filter the data to only include the selected columns
    df = (
        df
        .reindex(columns=[cname.split(" (")[0] for cname in columns])
        .dropna()
    )

    tabs["cluster"].write(
        f"Number of samples with all values: {df.shape[0]:,}"
    )

    # Filter to the samples which are being displayed here
    metadata = metadata.reindex(index=df.index)

    # Drop any columns which are entirely missing
    # Drop any columns which are invariant
    # Drop any columns which are entirely unique
    metadata = metadata.reindex(columns=[
        cname for cname, cvals in metadata.items()
        if (
            cvals.notnull().any() and 
            len(cvals.unique()) > 1 and 
            len(cvals.unique()) < metadata.shape[0]
        )
    ])

    # Let the user filter samples by metadata
    query_str = tabs["cluster"].text_input(
        "Filter samples by metadata",
        placeholder="e.g. 'cell_type == \"Neuron\"'"
    )
    if query_str:
        try:
            metadata = metadata.query(query_str)
            df = df.reindex(index=metadata.index)
            tabs["cluster"].write(
                f"Number of samples after filtering by metadata: {df.shape[0]:,}" # noqa
            )
        except Exception as e:
            tabs["cluster"].write(f"Error filtering metadata: {e}")
    else:
        tabs["cluster"].write("No filtering performed")

    tabs["cluster"].dataframe(metadata)

    return df, metadata


def format_col_labels(df: pd.DataFrame):

    # Count up the proportion of missing values in each column
    # and format a nice label
    col_labels = pd.DataFrame([
        {
            "column": cname,
            "missing": cvals.isnull().mean(),
            "label": f"{cname} ({round(100 * cvals.isnull().mean(), 1)}% missing)", # noqa
            "ix": ix
        }
        for ix, (cname, cvals) in enumerate(df.items())
    ])

    # By default, show the columns which are all tied for the least missing
    default_col_labels = col_labels.loc[
        col_labels["missing"] == col_labels["missing"].min()
    ]

    # Only include 10 by default (sorting by position in the input DataFrame)
    if len(default_col_labels) > 10:
        default_col_labels = default_col_labels.head(10)

    return col_labels["label"].tolist(), default_col_labels["label"].tolist()


@st.cache_data
def run_clustering(df: pd.DataFrame, **kwargs):
    clusters = (
        cluster.KMeans(**kwargs)
        .fit_predict(df.values)
    )
    return list(map(str, clusters))


@st.cache_data
def run_umap(df, **kwargs) -> pd.DataFrame:
    reducer = umap.UMAP(**kwargs)
    return pd.DataFrame(
        reducer.fit_transform(df),
        index=df.index,
        columns=["UMAP1", "UMAP2"]
    )


def setup_inputs():
    st.sidebar.file_uploader(
        label="Metadata (CSV)",
        type=["csv", "xlsx"],
        key="metadata_file",
        help="Left-most column should be the unique observation IDs."
    )

    for i in range(1, st.session_state["n_inputs"]+1):

        st.sidebar.text_input(
            "Measurement Type",
            value=f"Type {i}",
            key=f"mod{i}_type"
        )
        st.sidebar.file_uploader(
            f"Data Table {i}",
            key=f"mod{i}_file",
            type=["csv", "xlsx"],
            help="Left-most column should be the unique observation IDs."
        )


def setup_outputs():
    display = st.container()
    global tabs
    tabs = dict(zip(
        [
            "cluster",
            # "compare",
            "logs"
        ],
        display.tabs([
            "Cluster Analysis",
            # "Compare Modalities",
            "Data Ingest",
        ])
    ))


def read_files():
    for key in (
        ["metadata_file"] +
        [
            f"mod{i}_file"
            for i in range(1, st.session_state["n_inputs"] + 1)
        ]
    ):
        uploaded_file = st.session_state.get(key)
        if uploaded_file is not None:
            fn = uploaded_file.name
            parse_func = (
                pd.read_csv
                if fn.endswith(".csv")
                else pd.read_excel
            )
            df: pd.DataFrame = (
                parse_func(
                    uploaded_file,
                    index_col=0
                )
                .rename(
                    index=str
                )
            )
            df = df.reindex(columns=[
                cname for cname, cvals in df.items()
                if cvals.apply(lambda v: isinstance(v, (float, int))).all()
            ])
            st.session_state[key + '_data'] = df
            tabs["logs"].write(f"Read {key.replace('_', ' ')}: {fn} ({df.shape[0]:,} rows and {df.shape[1]:,} columns)") # noqa


def format_mudata():

    # Get the set of samples which have metadata defined
    metadata = st.session_state.get("metadata_file_data")
    if metadata is None:
        tabs["logs"].write("No metadata file provided")
        return
    obs = metadata.index

    # Make a dict with all of the data provided
    data = {}

    for i in range(1, st.session_state["n_inputs"]+1):
        modality = st.session_state[f"mod{i}_type"]
        df = st.session_state.get(f"mod{i}_file_data")
        if df is None:
            continue
        input_n = df.shape[0]
        overlap = set(df.index) & set(obs)
        overlap_n = len(overlap)
        tabs["logs"].write(f"Modality {modality} has {input_n:,} rows, {overlap_n:,} of which are in the metadata") # noqa
        if overlap_n > 0:
            data[modality] = AnnData(X=df.reindex(index=obs), obs=metadata)
            tabs["logs"].write(f"Adding modality {modality}: {data[modality].shape[1]:,} features") # noqa
        else:
            tabs["logs"].write(f"Skipping modality {modality} as there are no overlapping samples") # noqa

    if len(data) > 0:
        # Create MuData object
        st.session_state["mudata"] = MuData(data)


if __name__ == '__main__':
    main()
