import streamlit as st
import logging
import pandas as pd
from mudata import MuData
from anndata import AnnData
from mudata_explorer.app.mdata import set_mdata
from mudata_explorer.app.hash import get_dat_hash, set_mdata_hash
from mudata_explorer.helpers.io import hydrate_uns
from mudata_explorer.sdk import view

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def _read_mdata():
    # Ask the user to upload a CSV
    st.write("Please upload the CSV file with the viral titer data.")
    uploaded_file = st.file_uploader("Choose a file")
    if uploaded_file is None:
        return
    
    # Read the CSV
    df = pd.read_csv(uploaded_file)

    # Make a dict with the wide dataframes
    adata_dict = {
        label: adata
        for label, adata in _parse_df(df)
    }

    # Ask the user which measurements to include
    to_keep = st.multiselect(
        "Select the measurements to include:",
        list(adata_dict.keys()),
        default=[]
    )

    if not to_keep:
        return

    mdata = MuData({
        label: adata_dict[label]
        for label in to_keep
    })

    mdata.obs = list(adata_dict.values())[0].obs

    # Filter the data based on visit days
    mdata = _filter_mdata(mdata, "visit days", "visitno")
    if mdata is None:
        return

    # Filter the data based on treatment arms
    mdata = _filter_mdata(mdata, "treatment arm", "TRT01P")
    if mdata is None:
        return
    
    return mdata


def _filter_mdata(mdata: MuData, label: str, column: str):
    selections = st.multiselect(
        f"Select the {label} to include:",
        mdata.obs[column].unique(),
        default=[]
    )

    if not selections:
        return

    ix = mdata.obs[column].isin(selections)

    return mdata[ix]


def _add_figures(mdata: MuData):

    view.markdown(
        mdata,
        text_value="## Viral Titer Analysis (HVTN Format)",
    )

    for variable_name in mdata.mod.keys():

        # Get the values to show
        df = mdata.mod[variable_name].to_df()
        values_to_show = [
            cname
            for cname, cvals in df.items()
            if cvals.std() > 0
        ]

        for cname in values_to_show:

            _make_boxplot(mdata, variable_name, cname, x="visitchar", color="TRT01P")
            _make_boxplot(mdata, variable_name, cname, color="visitchar", x="TRT01P")

def _make_boxplot(mdata: MuData, variable_name: str, cname: str, x="visitchar", color="TRT01P"):

    view.plotly_box(
        mdata,
        **{
            "data": {
                "sidebar": False,
                "axis": {
                    "value": 0,
                    "sidebar": False
                },
                "transforms": {
                    "value": [],
                    "sidebar": False
                },
                "columns": {
                    "x": {
                        "sidebar": False,
                        "table": {
                            "value": "Observation Metadata",
                            "sidebar": False
                        },
                        "cname": {
                            "value": x,
                            "sidebar": False
                        },
                        "label": {
                            "value": x,
                            "sidebar": False
                        },
                        "scale": {
                            "value": None,
                            "sidebar": False
                        },
                        "colorscale": False,
                        "is_categorical": {
                            "value": False,
                            "sidebar": False
                        }
                    },
                    "y": {
                        "sidebar": False,
                        "table": {
                            "value": f"{variable_name}.data",
                            "sidebar": False
                        },
                        "cname": {
                            "value": cname,
                            "sidebar": False
                        },
                        "label": {
                            "value": cname,
                            "sidebar": False
                        },
                        "scale": {
                            "value": None,
                            "sidebar": False
                        },
                        "colorscale": False,
                        "is_categorical": {
                            "value": False,
                            "sidebar": False
                        }
                    },
                    "color": {
                        "enabled": {
                            "value": True,
                            "sidebar": False
                        },
                        "sidebar": False,
                        "table": {
                            "value": "Observation Metadata",
                            "sidebar": False
                        },
                        "cname": {
                            "value": color,
                            "sidebar": False
                        },
                        "label": {
                            "value": color,
                            "sidebar": False
                        },
                        "scale": {
                            "value": "D3",
                            "sidebar": False
                        },
                        "colorscale": True,
                        "is_categorical": {
                            "value": True,
                            "sidebar": False
                        }
                    }
                },
                "filter_rows": {
                    "sidebar": False,
                    "type": {
                        "value": None,
                        "sidebar": False
                    },
                    "tables": {
                        "value": [],
                        "sidebar": False
                    },
                    "cname": {
                        "value": None,
                        "sidebar": False
                    },
                    "expr": {
                        "value": None,
                        "sidebar": False
                    },
                    "value_enum": {
                        "value": None,
                        "sidebar": False
                    },
                    "value_str": {
                        "value": None,
                        "sidebar": False
                    }
                }
            },
            "scale_options": {
                "log_y": {
                    "value": None,
                    "sidebar": True
                }
            },
            "formatting": {
                "title": {
                    "value": f"{variable_name} -- {cname}",
                    "sidebar": True
                },
                "legend": {
                    "value": None,
                    "sidebar": False
                }
            },
            "statistics": {
                "compare_groups": {
                    "value": "Disabled",
                    "sidebar": True
                }
            }
        }
    )




def _parse_df(df):

    # Column Categories:
    mod_columns = [
        "titer_crit"
    ]

    val_columns = [
        "titer_num",
        "neutralization_pct_init_dilution"
    ]

    obs_columns = [
        "ptid",
        "visitno",
        "incubation",
        "experimenter",
        "TRT01P",
        "panel",
        "visitchar",
        "visitmon",
        "TRT01PG",
        "TRT01P_P",
        "platenum"
    ]
    
    var_columns = [
        "isolate"
    ]
    
    ign_columns = [
        "virusdilution",
        "poscrit",
        "out_win",
        "not_exp",
        "r_titer",
        "tier",
        "virusid",
    ]

    # Drop the ignore columns
    df = df.drop(columns=ign_columns)

    # Iterate over every value column
    for col in val_columns:
        logger.info(f"Processing column {col}")
        # Iterate over each combination of the mod columns
        for group_ids, group_df in df.groupby(mod_columns):
            titer = group_ids[0]
            logger.info(f"Processing group {group_ids}")
            
            # Make a wide table
            wide_df = (
                group_df.pivot(
                    index=obs_columns,
                    columns=var_columns,
                    values=col
                )
                .rename(
                    columns=lambda x: (
                        x
                        .replace("\\", "_")
                        .replace("__", "_")
                    )
                )
            )
            obs = wide_df.index.to_frame(index=False)

            # Make an AnnData object
            adata = AnnData(
                X=wide_df.reset_index(drop=True),
                obs=obs
            )

            label = f"{col} {titer}"

            yield label, adata


def _save_data(mdata: MuData):
    # Get the hash of the data
    _, hash, _ = get_dat_hash(mdata)

    hydrate_uns(mdata)
    set_mdata(mdata, full=True, id="main")
    set_mdata_hash(hash, id="main")
    st.switch_page("pages/views.py")


def main():
    st.write("## Viral Titer Analysis (HVTN Format)")

    # Read in the data
    mdata = _read_mdata()
    if mdata is None:
        return
    
    if not st.button("Build Figures"):
        return

    # Add the figures
    _add_figures(mdata)

    # Save the data
    _save_data(mdata)

    # Navigate to the figures
    st.switch_page("pages/view_all.py")


if __name__ == "__main__":
    main()
