from functools import lru_cache
from io import BytesIO
import anndata as ad
from cirro import DataPortalDataset
from cirro.sdk.file import DataPortalFile # noqa
from cirro.models.file import File
import streamlit as st
from streamlit.delta_generator import DeltaGenerator
from tempfile import TemporaryDirectory
from time import sleep
from typing import List, Optional
import re
import pandas as pd
from mudata import MuData
from mudata_explorer.sdk import process, view


def clear_cirro_client():
    # Wait for a second
    sleep(1)
    # Clear the client
    del st.session_state["Cirro"]
    # Let the user log in again
    st.rerun()


def list_files(
    dataset: DataPortalDataset,
    pattern: str = None
) -> List[DataPortalFile]:
    """List the files in a dataset."""

    # Try to get the list of files in the dataset
    try:
        files = dataset.list_files()
    # If there is an error
    except Exception as e:
        # Report it to the user
        st.exception(e)
        clear_cirro_client()

    # If a pattern was provided
    if pattern is not None:
        # Filter the files
        files = [f for f in files if re.match(pattern, f.name[5:])]

    return files


@lru_cache
def process_name(process_id: str) -> str:
    """Return the readable name of a process defined in Cirro."""

    # Return the name attribute of the process with the given ID
    return (
        st.session_state["Cirro"]
        .get_process_by_id(process_id)
        .name
    )


def find_file_by_extension(
    files: List[DataPortalFile],
    prefix="data/",
    suffix="",
    selectbox_label="Use file:",
    none_found_msg="No file found"
):
    # Filter to only the relative abundance tables
    files = [
        f
        for f in files
        if f.name.startswith(prefix) and f.name.endswith(suffix)
    ]

    # If there are no relative abundance tables
    if len(files) == 0:
        st.error(none_found_msg)
        return
    # If there is more than one relative abundance table
    elif len(files) > 1:
        # Pick the file from the dataset
        selected_abund = st.selectbox(
            selectbox_label,
            [f.name[len(prefix):-len(suffix)] for f in files],
            index=None,
            placeholder="< select a file >"
        )
        if selected_abund == "< select a file >":
            return
        for file in files:
            if file.name == f"{prefix}{selected_abund}{suffix}":
                return file
    else:
        return files[0]


def select_file(
    container: DeltaGenerator,
    dataset: DataPortalDataset,
    pattern: str = None,
    files: Optional[List[DataPortalFile]] = None
) -> DataPortalFile:

    # If no files were provided
    if files is None:
        # Try to get the list of files in the dataset
        files = list_files(container, dataset)
    else:
        assert all(isinstance(f, DataPortalFile) for f in files)

    # Give the user the option to select a file
    # (strip off the 'data/' prefix)
    file = container.selectbox(
        "File",
        ["< select a file >"] + [
            f.name[5:] for f in files
            if pattern is None or re.match(pattern, f.name[5:])
        ],
        key="select-file"
    )

    if file == "< select a file >":
        return

    # Return the file object
    return next(f for f in files if f.name[5:] == file)


def sample_metadata(dataset: DataPortalDataset) -> pd.DataFrame:
    """Return the samplesheet saved with a dataset."""
    access_context = dataset.list_files()[0]._file.access_context
    return pd.read_csv(
        BytesIO(
            dataset._client._file_service.get_file(
                File(
                    relative_path="config/samplesheet.csv",
                    size=0,
                    access_context=access_context
                )
            )
        ),
        index_col=0
    )


@st.cache_data(hash_funcs={DataPortalFile: lambda f: f.absolute_path})
def read_h5ad(h5ad_file: DataPortalFile) -> ad.AnnData:
    # Download the file to a temporary directory
    with TemporaryDirectory() as tmp:
        h5ad_file.download(tmp)
        filename = f"{tmp}/{h5ad_file.name}"
        return ad.read_h5ad(filename)


def _guess_if_categorical(vals: pd.Series) -> bool:
    """
    Try to infer from a collection of values whether it is
    categorical or continuous in nature.
    True if categorical else False
    """
    vals = vals.dropna()
    if vals.apply(lambda v: isinstance(v, str)).any():
        return True

    vc = vals.value_counts()

    # If more than 25% of the values are singletons
    if (vc == 1).sum() > (0.25 * vals.shape[0]):
        return False

    return True


def ask_if_categorical(kw: str, vals: pd.Series) -> bool:
    guess_is_categorical = _guess_if_categorical(vals)
    guess_str = 'categorical' if guess_is_categorical else 'continuous'

    st.write(f"""**Categorical / Continuous**

The variable '{kw}' appears to be {guess_str}.
The user may indicate whether it should be treated in this way for display.""")
    return st.selectbox(
        f"Variable Type ('{kw}')",
        ["Continuous", "Categorical"],
        index=int(guess_is_categorical)
    ) == "Categorical"


def summarize_features(
    mdata: MuData,
    mod="mod"
):
    process.summary_stats(
        mdata,
        table_data_tables_value=f"{mod}.data",
        table_data_axis_value=1
    )


def shannon_diversity(
    mdata: MuData,
    mod="mod",
    dest_key="shannon"
):
    process.shannon_diversity(
        mdata,
        outputs_dest_key_value=dest_key,
        table_data_axis_value=0,
        table_data_tables_value=[f"{mod}.data"]
    )


def run_umap(
    mdata: MuData,
    mod="mod",
    metric="braycurtis",
    n_neighbors=15,
    min_dist=0.5,
    dest_key="umap",
    n_components=2
):
    process.umap(
        mdata,
        umap_params_metric_value=metric,
        umap_params_n_neighbors_value=n_neighbors,
        umap_params_min_dist_value=min_dist,
        umap_params_n_components_value=n_components,
        outputs_dest_key_value=dest_key,
        table_data_axis_value=0,
        table_data_tables_value=[f"{mod}.data"]
    )


def run_pca(
    mdata: MuData,
    mod="mod",
    dest_key="pca"
):
    process.pca(
        mdata,
        outputs_dest_key_value=dest_key,
        table_data_axis_value=0,
        table_data_tables_value=[f"{mod}.data"]
    )


def run_leiden(
    mdata: MuData,
    mod="mod",
    metric="braycurtis",
    resolution=1.0,
    n_neighbors=15,
    dest_key="leiden",
):
    process.leiden(
        mdata,
        clustering_metric_value=metric,
        clustering_n_neighbors_value=n_neighbors,
        clustering_resolution_value=resolution,
        outputs_dest_key_value=dest_key,
        table_data_axis_value=0,
        table_data_tables_value=[f"{mod}.data"]
    )
    assert dest_key in mdata.obs.columns, mdata


def run_kmeans(
    mdata: MuData,
    mod="mod",
    k=15,
    transforms=[],
    dest_key="kmeans"
):
    # Automatically check if max_k is greater than the number of samples
    n_samples = mdata.shape[0]
    if n_samples < 10:
        max_k = n_samples
    else:
        max_k = max(k, 10)
    if n_samples < k:
        k = n_samples

    process.kmeans(
        mdata,
        clustering_k_value=k,
        outputs_dest_key_value=dest_key,
        table_data_axis_value=0,
        table_data_tables_value=[f"{mod}.data"],
        table_data_transforms_value=transforms,
        clustering_max_k_value=max_k
    )


def run_kruskal(
    mdata: MuData,
    table: str,
    dest_key: str,
    grouping_cname: str,
    grouping_table: str
):
    process.kruskal(
        mdata,
        table_data_axis_value=0,
        table_data_tables_value=[table],
        table_data_transforms_value=[],
        outputs_dest_key_value=dest_key,
        table_grouping_axis_value=0,
        table_grouping_columns_grouping_cname_value=grouping_cname,
        table_grouping_columns_grouping_label_value=grouping_cname,
        table_grouping_columns_grouping_table_value=[grouping_table]
    )


def run_deseq2(
    mdata: MuData,
    table: str,
    dest_key: str,
    grouping_cname: str,
    grouping_table: str,
    ref_level: str,
    comp_level: str,
    alpha=0.05,
    independent_filter=True
):
    process.deseq2(
        mdata,
        table_data_axis_value=0,
        table_data_tables_value=[table],
        table_data_transforms_value=[],
        outputs_dest_key_value=dest_key,
        table_grouping_axis_value=0,
        table_grouping_columns_grouping_cname_value=grouping_cname,
        table_grouping_columns_grouping_label_value=grouping_cname,
        table_grouping_columns_grouping_table_value=[grouping_table],
        comparison_ref_level_value=ref_level,
        comparison_comp_level_value=comp_level,
        comparison_alpha_value=alpha,
        comparison_independent_filter_value=independent_filter
    )


def run_spearman(
    mdata: MuData,
    table: str,
    dest_key: str,
    comparitor_cname: str,
    comparitor_table: str
):
    process.spearman(
        mdata,
        outputs_dest_key_value=dest_key,
        table_data_axis_value=0,
        table_data_tables_value=[table],
        table_data_transforms_value=[],
        table_comparitor_axis_value=0,
        table_comparitor_columns_comparitor_cname_value=comparitor_cname,
        table_comparitor_columns_comparitor_label_value=comparitor_cname,
        table_comparitor_columns_comparitor_table_value=[comparitor_table]
    )


def add_scatter(
    mdata: MuData,
    title: str,
    legend: str,
    x: str,
    xlabel: str,
    y: str,
    ylabel: str,
    table: str,
    axis: int,
    color_table: str,
    cname: str,
    label: str,
    is_categorical: bool,
    scale: str,
    xtable=None,
    ytable=None,
    log_x=False,
    log_y=False
):

    view.plotly_scatter(
        mdata,
        formatting_title_value=title,
        formatting_title_sidebar=True,
        data_axis_value=axis,
        data_columns_color_cname_value=cname,
        data_columns_color_is_categorical_value=is_categorical,
        data_columns_color_enabled_value=color_table is not None,
        data_columns_color_label_value=label,
        data_columns_color_scale_value=scale,
        data_columns_color_table_value=[color_table],
        data_columns_size_enabled_value=False,
        data_columns_x_cname_value=x,
        data_columns_x_cname_sidebar=True,
        data_columns_x_label_value=xlabel,
        data_columns_x_table_value=[table if xtable is None else xtable],
        data_columns_y_cname_value=y,
        data_columns_y_cname_sidebar=True,
        data_columns_y_label_value=ylabel,
        data_columns_y_table_value=[table if ytable is None else ytable],
        formatting_legend_value=legend,
        formatting_legend_sidebar=True,
        scale_options_log_x_value=log_x,
        scale_options_log_y_value=log_y,
        scale_options_log_x_sidebar=False,
        scale_options_log_y_sidebar=False,
        formatting_opacity_sidebar=False
    )


def add_stacked_bars(
    mdata: MuData,
    title: str,
    yaxis_title: str,
    table: str,
    category_cname: str,
    category_label: str,
    category_table="Observation Metadata",
    sort_rows_by="Values",
    features=None,
    legend: Optional[str]=None
):
    # If a set of features was specified
    if features is not None:
        # Filter the columns to only include those features
        filter_kwargs = dict(
            table_data_filter_cols_tables_value=[table],
            table_data_filter_cols_type_value="index",
            table_data_filter_cols_expr_value="in",
            table_data_filter_cols_value_enum_value=features,
            table_data_filter_cols_value_enum_sidebar=True
        )
    else:
        filter_kwargs = dict()

    view.plotly_stacked_bars(
        mdata,
        formatting_title_value=title,
        formatting_yaxis_title_value=yaxis_title,
        formatting_sort_rows_by_value=sort_rows_by,
        table_category_axis_value=0,
        table_category_columns_category_cname_value=category_cname,
        table_category_columns_category_is_categorical_value=True,
        table_category_columns_category_label_value=category_label,
        table_category_columns_category_table_value=[category_table],
        table_category_enabled_value=True,
        table_data_tables_value=[table],
        formatting_legend_value=legend,
        **filter_kwargs
    )


def add_table(
    mdata: MuData,
    title: str,
    legend: str,
    table: str,
    sort_by: str,
    axis: int,
):
    view.table(
        mdata,
        options_sort_columns_sort_by_cname_value=sort_by,
        options_sort_columns_sort_by_label_value=sort_by,
        options_sort_columns_sort_by_table_value=[table],
        options_sort_axis_value=axis,
        data_table_axis_value=axis,
        data_table_tables_value=[table],
        display_options_title_value=title,
        display_options_legend_value=legend
    )


def add_boxplot(
    mdata: MuData,
    title: str,
    legend: str,
    x_table: str,
    x_cname: str,
    x_label: str,
    y_table: str,
    y_cname: str,
    y_label: str,
):
    view.plotly_box(
        mdata,
        formatting_title_value=title,
        statistics_compare_groups_value="Kruskal-Wallis",
        data_columns_x_table_value=[x_table],
        data_columns_x_cname_value=x_cname,
        data_columns_x_cname_sidebar=True,
        data_columns_x_label_value=x_label,
        data_columns_y_table_value=[y_table],
        data_columns_y_cname_value=y_cname,
        data_columns_y_cname_sidebar=True,
        data_columns_y_label_value=y_label,
        data_columns_color_enabled_value=False,
        formatting_legend_value=legend,
        formatting_legend_sidebar=True,
        scale_options_log_y_sidebar=False,
    )


def add_boxplot_multi(
    mdata: MuData,
    mod="mod",
    cols_query={},
    rows_query={},
    var_label="Feature",
    val_label="Value",
    title="",
    legend=""
):
    view.plotly_box_multiple(
        mdata,
        table_category_enabled_value=False,
        table_data_tables_value=[f"{mod}.data"],
        table_data_filter_cols=cols_query,
        table_data_filter_rows=rows_query,
        variable_options_sort_by_value="Median",
        variable_options_axis_value="X-Axis",
        display_options_val_label_value=val_label,
        display_options_var_label_value=var_label,
        display_options_title_value=title,
        display_options_legend_value=legend
    )


def add_category_count(
    mdata: MuData,
    title: str,
    legend: str,
    x_table: str,
    x_cname: str,
    x_label: str,
    color_table: str,
    color_cname: str,
    color_label: str
):
    view.plotly_category_count(
        mdata,
        data_columns_x_table_value=[x_table],
        data_columns_x_cname_value=x_cname,
        data_columns_x_label_value=x_label,
        data_columns_color_table_value=[color_table],
        data_columns_color_cname_value=color_cname,
        data_columns_color_label_value=color_label,
        barmode_value="group",
        formatting_title_value=title,
        annotation_options_chisquare_value=True,
        formatting_legend_value=legend
    )


def format_float(f: float) -> str:
    if f > 0.01:
        return f"{f:.2f}"
    elif f > 0.001:
        return f"{f:.3f}"
    elif f > 0.0001:
        return f"{f:.4f}"
    else:
        return f"{f:.2e}"
