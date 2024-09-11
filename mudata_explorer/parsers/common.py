from anndata import AnnData
from dataclasses import dataclass
from mudata_explorer.sdk import io
from mudata_explorer.helpers.cirro_readers import util
from mudata import MuData
import pandas as pd
import streamlit as st
from typing import Optional


@dataclass
class CompositionalParams:
    """Configuration for the analysis of compositional data."""
    # Name of the dataset, to be used in the title
    dataset_name: str = "Compositional Analysis"
    # Column name to use for comparison
    compare_by: Optional[str] = None
    # Label to use when describing those comparison groups
    label: Optional[str] = None
    # Whether the comparison variable is categorical
    is_categorical: Optional[str] = None
    # The comparison may be categorical, optionally
    force_categorical: bool = False
    # Number of features to use in the analysis and plots
    n_top_features: int = 20

    @classmethod
    def get_params(cls, mdata: MuData):

        params = dict()

        # If there is sample metadata
        if mdata.obs.shape[1] > 0:

            # Let the user pick a column to compare samples by
            with st.container(border=1):
                st.write("""**Compare Samples By:**

The user may optionally select a column from the sample metadata to
compare samples by.
The comparison will depend on whether the variable is continuous
(like age or weight) or categorical (discrete groups).
    """)
                compare_by = st.selectbox(
                    "Sample Metadata",
                    ["< select a column >"] + list(mdata.obs.columns),
                    index=0
                )

                if compare_by != "< select a column >":
                    params["compare_by"] = compare_by

                    # Get a readable label for this column
                    st.write("**Variable Label**")
                    params["label"] = st.text_input(
                        f"Label used in displays for '{compare_by}':",
                        compare_by.replace("_", " ")
                    )

                    if not cls.force_categorical:
                        params["is_categorical"] = util.ask_if_categorical(
                            params["compare_by"],
                            mdata.obs[params["compare_by"]],
                        )

                    if cls.force_categorical or params["is_categorical"]:
                        if hasattr(cls, "ref_group"):
                            params["ref_group"] = st.selectbox(
                                "Reference Group",
                                list(mdata.obs[params["compare_by"]].unique()),
                                index=None
                            )
                        if hasattr(cls, "comp_group"):
                            params["comp_group"] = st.selectbox(
                                "Comparison Group",
                                [
                                    val for val in mdata.obs[params["compare_by"]].unique()
                                    if val != params["ref_group"]
                                ],
                                index=None
                            )
                        if hasattr(cls, "deseq2_alpha"):
                            params["deseq2_alpha"] = st.number_input(
                                "DESeq2 Alpha",
                                min_value=0.01,
                                max_value=0.1,
                                value=cls.deseq2_alpha,
                                step=0.01
                            )

        if hasattr(cls, "n_top_features"):
            with st.container(border=1):
                st.write("""**Top Features**

An extended set of figures will be displayed for the figures which are
the most strongly associated with the variable of interest.
The user may decide how many such variables to include in the display.
        """)
                params["n_top_features"] = st.number_input(
                    "Number of features to display (with the lowest p-values)",
                    min_value=1,
                    max_value=100,
                    value=10,
                    step=1
                )

        if hasattr(cls, "leiden_res"):
            with st.container(border=1):
                st.write("""**Sample Clustering**

The samples will be clustered into groups to identify
any subpopulations that share similar patterns of features.
                            """)
                params["leiden_res"] = st.number_input(
                    "Leiden Resolution - Higher values lead to more clusters",
                    min_value=0.1,
                    max_value=2.0,
                    value=1.0,
                    step=0.1
                )

        if hasattr(cls, "k"):
            with st.container(border=1):
                st.write("""**Sample Clustering**

The samples will be clustered into groups to identify
any subpopulations that share similar patterns of features.
                            """)
                params["k"] = st.number_input(
                    "K-Means Clustering - Number of clusters",
                    min_value=1,
                    max_value=100,
                    value=int(cls.k)
                )

        if hasattr(cls, "plot_genes"):
            with st.container(border=1):
                st.write("""**Plot Specific Genes**

Provide a list of any specific genes which should be plotted
across the sample groups.
                            """)
                params["plot_genes"] = ",".join(st.multiselect(
                    "Genes to plot",
                    list(mdata.mod["expression"].var.index),
                    [
                        val
                        for val in cls.plot_genes.split(",")
                        if val in mdata.mod["expression"].var.index
                    ]
                ))

        params["dataset_name"] = st.text_input(
            "Dataset Name",
            cls.dataset_name
        )

        return cls(**params)


def top_features(
    mdata: MuData,
    n_top_features: int,
    mod="abund",
    varm_key="summary_stats",
    cname="mean",
    ascending=False
):
    return list(
        mdata
        .mod[mod]
        .varm[varm_key]
        .sort_values(
            by=cname,
            ascending=ascending
        )
        .head(n_top_features)
        .index.values
    )


def top_features_filter_cols(
    mdata: MuData,
    n_top_features: int,
    mod="abund",
    varm_key="summary_stats",
    cname="mean",
    ascending=False
): 
    features = top_features(
        mdata,
        n_top_features,
        varm_key=varm_key,
        cname=cname,
        ascending=ascending,
        mod=mod
    )

    return dict(
        tables_value=[f"{mod}.data"],
        type_value="index",
        expr_value="in",
        value_enum_value=features,
        value_enum_sidebar=True
    )


def parse_adata(
    adata: AnnData,
    groupby_var=False,
    sum_to_one=True,
    kw="abund"
) -> MuData:

    # Filter out any samples with 0 counts
    adata = adata[adata.to_df().sum(axis=1) > 0]
    assert adata.shape[0] > 0, "No samples with non-zero counts."

    if sum_to_one:
        # Make sure that every row sums to 1
        adata.X = adata.to_df().apply(lambda r: r / r.sum(), axis=1).values

    if groupby_var:
        # Pick a taxonomic level to collapse by
        with st.container(border=1):
            adata = _collapse_by_taxon(adata)

    # Get the sample metadata
    with st.container(border=1):
        adata = _read_sample_meta(adata)

    # Optionally filter samples by sample metadata
    with st.container(border=1):
        adata = _filter_samples(adata)

    # Coerce objects to str in the metadata table
    adata.obs = adata.obs.astype({
        cname: "str"
        for cname, dtype in adata.obs.dtypes.items()
        if dtype == "object"
    })

    # Make a MuData object
    mdata = io.build_mdata(
        {
            kw: AnnData(
                X=adata.to_df(),
                var=adata.var
            )
        },
        obs=adata.obs
    )

    return mdata


def _collapse_by_taxon(adata: AnnData) -> AnnData:
    # The taxonomy is in the var table. However, it is not labelled.
    filt_levels = list(adata.var.columns.values)
    assert len(filt_levels) > 0

    # Let the user select which taxon to collapse by
    level = st.selectbox(
        "Label organisms to the level of:",
        filt_levels,
        index=max(0, len(filt_levels)-2)
    )

    # If the ASV isn't classified at that level, fill in with the
    # next highest level
    adata.var["_group_level"] = adata.var.apply(
        lambda r: _pick_level(r, filt_levels.index(level)),
        axis=1
    )

    # Sum up the abundances for each taxon
    abund = (
        adata
        .to_df()
        .T
        .groupby(adata.var["_group_level"])
        .sum()
        .sort_index()
        .rename_axis(index=level)
        .T
    )

    # Show the user the most abundant taxa at this level
    st.dataframe(
        abund
        .mean()
        .sort_values(ascending=False)
        .reset_index()
        .rename(columns={0: "Average Relative Abundance"}),
        hide_index=True
    )

    # Get the taxonomic assignment for each taxon
    tax = (
        adata.var
        .rename(
            columns={
                f"taxonomy_{ix}": label
                for ix, label in enumerate(filt_levels)
            }
        )
        .reindex(columns=filt_levels[:filt_levels.index(level)] + ["_group_level"])
        .groupby("_group_level")
        .head(1)
        .set_index("_group_level")
        .sort_index()
        .rename_axis(index=level)
        .reindex(index=abund.columns.values)
    )

    # Return an AnnData object
    return AnnData(X=abund, var=tax, obs=adata.obs)


def _pick_level(row: pd.Series, pick_ix: int):

    # Walk from each of the levels back to the first
    for ix in range(pick_ix, -1, -1):

        # If there is a value at this level
        if not pd.isnull(row.iloc[ix]):
            return row.iloc[ix]

    # If none of the levels have values
    return "Unclassified"


def _read_sample_meta(adata: AnnData) -> AnnData:

    st.write("""
**Sample Metadata**

The table below shows the metadata associated with each sample in the dataset.
To modify or augment this information, download as a CSV and upload a new copy.
""")

    # Let the user view the metadata
    st.write(adata.obs.reset_index())

    # Provide a button to download as CSV
    st.download_button(
        label="Download sample metadata",
        data=adata.obs.to_csv(),
        file_name="sample_metadata.csv",
        mime="text/csv"
    )

    # Let the user upload a new samplesheet
    uploaded_file = st.file_uploader("Upload a new spreadsheet (in CSV or XLSX format)")

    if uploaded_file:
        if uploaded_file.name.endswith(".csv"):
            new_meta = pd.read_csv(uploaded_file, index_col=0)
        elif uploaded_file.name.endswith(".xlsx"):
            new_meta = pd.read_excel(uploaded_file, index_col=0)
        else:
            st.error("Please upload a CSV or XLSX file.")
            return

        # Subset to only those rows which are shared
        n_intersection = new_meta.index.intersection(adata.obs.index).shape[0]
        st.write(f"The uploaded file has {n_intersection:,} valid sample labels.")
        if n_intersection > 0:

            # Drop any columns which are entirely missing
            new_meta = new_meta.drop(
                columns=[
                    cname
                    for cname in new_meta.columns
                    if new_meta[cname].isnull().all()
                ]
            )

            # Show the user the updated metadata
            new_meta = new_meta.reindex(adata.obs.index)
            st.dataframe(new_meta)
            adata.obs = new_meta

    return adata


def _filter_samples(adata: AnnData):

    st.markdown("""**Filter Samples**

Optionally filter the set of samples which are included in the analysis
""")
    
    # Let the user select the approach for filtering samples
    filter_method = st.selectbox(
        "Filter Samples",
        [
            "Query String",
            "Select Samples to Include",
            "Select Samples to Exclude"
        ],
        index=None,
        placeholder="Select a method"
    )
    if filter_method == "Query String":
        return _filter_samples_query(adata)
    elif filter_method == "Select Samples to Include":
        return _filter_samples_select(adata, include=True)
    elif filter_method == "Select Samples to Exclude":
        return _filter_samples_select(adata, include=False)
    else:
        return adata


def _filter_samples_query(adata: AnnData) -> AnnData:
    ntot = adata.shape[0]
    
    query_string = st.text_input(
        "Query String",
        placeholder="e.g. group == 'A' or timepoint == 1"
    )
    if query_string is not None and len(query_string) > 0:
        try:
            filtered = adata[
                adata.obs.query(query_string).index.values
            ]
        except Exception as e:
            st.exception(e)

        if filtered.shape[0] > 0:

            st.markdown(f"Including {filtered.shape[0]:,} / {ntot:,} samples")
            return filtered
        
        else:
            st.write("No samples matching the query")

    else:
        st.markdown(f"Including all {ntot:,} samples")

    return adata


def _filter_samples_select(adata: AnnData, include: bool) -> AnnData:
    ntot = adata.shape[0]

    # Let the user select samples to include or exclude
    selected = st.multiselect(
        f"Select Samples to {'Include' if include else 'Exclude'}",
        adata.obs.index.values,
        adata.obs.index.values if include else []
    )

    if len(selected) > 0 and len(selected) < ntot:

        if include:
            filtered = adata[selected]
            st.markdown(f"Including {filtered.shape[0]:,} / {ntot:,} samples")
        else:
            filtered = adata[~adata.obs.index.isin(selected)]
            st.markdown(f"Including {filtered.shape[0]:,} / {ntot:,} samples")

        return filtered

    else:
        st.markdown(f"Including all {ntot:,} samples")
        return adata