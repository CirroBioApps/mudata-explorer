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
from muon import MuData
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
        st.error(f"Error: {e}")
        clear_cirro_client()

    # If a pattern was provided
    if pattern is not None:
        # Filter the files
        files = [f for f in files if re.match(pattern, f.name[5:])]

    return files


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
        table_data_tables=f"{mod}.data",
        table_data_axis=1
    )


def run_umap(
    mdata: MuData,
    mod="mod",
    metric="braycurtis",
    n_neighbors=15,
    min_dist=0.5,
    dest_key="umap",
):
    process.umap(
        mdata,
        **{
            "umap_params": {
                "metric": metric,
                "n_neighbors": n_neighbors,
                "min_dist": min_dist,
                "n_components": 2
            },
            "outputs": {
                "dest_key": dest_key
            },
            "table": {
                "data": {
                    "axis": 0,
                    "tables": [
                        f"{mod}.data"
                    ]
                }
            }
        }
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
        **{
            "clustering": {
                "metric": metric,
                "n_neighbors": n_neighbors,
                "resolution": resolution
            },
            "outputs": {
                "dest_key": dest_key
            },
            "table": {
                "data": {
                    "axis": 0,
                    "tables": [
                        f"{mod}.data"
                    ]
                }
            }
        }
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
        **{
            "outputs": {
                "dest_key": dest_key
            },
            "table": {
                "data": {
                    "axis": 0,
                    "tables": [table],
                    "transforms": []
                },
                "grouping": {
                    "axis": 0,
                    "grouping": {
                        "cname": grouping_cname,
                        "label": grouping_cname,
                        "table": [grouping_table]
                    }
                }
            }
        }
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
):

    view.plotly_scatter(
        mdata,
        formatting_title=title,
        **{
            "data": {
                "axis": axis,
                "color": {
                    "cname": cname,
                    "is_categorical": is_categorical,
                    "enabled": True,
                    "label": label,
                    "scale": scale,
                    "table": [color_table]
                },
                "size": {
                    "enabled": False
                },
                "x": {
                    "cname": x,
                    "label": xlabel,
                    "table": [table]
                },
                "y": {
                    "cname": y,
                    "label": ylabel,
                    "table": [table]
                }
            }
        }
    )

    view.markdown(mdata, text=legend)


def add_table(
    mdata: MuData,
    title: str,
    legend: str,
    table: str,
    sort_by: str,
    axis: int,
):
    if len(title) > 0:
        view.markdown(mdata, text=title)

    view.table(
        mdata,
        **{
            "options": {
                "sort": {
                    "sort_by": {
                        "table": [table],
                        "cname": sort_by,
                        "label": sort_by
                    },
                    "axis": axis
                }
            },
            "data": {
                "table": {
                    "axis": axis,
                    "tables": [table]
                }
            }
        }
    )

    if len(legend) > 0:
        view.markdown(mdata, text=legend)


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
        formatting_title=title,
        statistics_compare_groups="Kruskal-Wallis",
        **{
            "data": {
                "x": {
                    "table": [
                        x_table
                    ],
                    "cname": x_cname,
                    "label": x_label
                },
                "y": {
                    "table": [
                        y_table
                    ],
                    "cname": y_cname,
                    "label": y_label
                },
                "color": {
                    "enabled": False
                }
            }
        }
    )

    if len(legend) > 0:
        view.markdown(mdata, text=legend)


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
        table_category_enabled=False,
        **{
            "table": {
                "data": {
                    "tables": [
                        f"{mod}.data"
                    ],
                    "cols_query": {
                        "query": cols_query
                    },
                    "rows_query": {
                        "query": rows_query
                    }
                }
            },
            "variable_options": {
                "sort_by": "Median",
                "axis": "X-Axis"
            },
            "display_options": {
                "var_label": var_label,
                "val_label": val_label,
                "title": title
            }
        }
    )

    if len(legend) > 0:
        view.markdown(mdata, text=legend)


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
        **{
            "data": {
                "x": {
                    "table": [
                        x_table
                    ],
                    "cname": x_cname,
                    "label": x_label
                },
                "color": {
                    "table": [
                        color_table
                    ],
                    "cname": color_cname,
                    "label": color_label
                }
            },
            "barmode": "group",
            "formatting_title": title,
            "annotation_options": {
                "chisquare": True
            }
        }
    )

    if len(legend) > 0:
        view.markdown(mdata, text=legend)


def format_float(f: float) -> str:
    if f > 0.01:
        return f"{f:.2f}"
    elif f > 0.001:
        return f"{f:.3f}"
    elif f > 0.0001:
        return f"{f:.4f}"
    else:
        return f"{f:.2e}"
