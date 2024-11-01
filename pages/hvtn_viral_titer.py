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

    mdata = MuData(adata_dict)
    mdata.obs = list(adata_dict.values())[0].obs

    return mdata


def _add_figures(mdata: MuData):

    view.markdown(
        mdata,
        text_value="## Viral Titer Analysis (HVTN Format)",
    )

    for variable_name in mdata.mod_names:

        # Get the values to show
        df = mdata.mod[variable_name].to_df()
        values_to_show = [
            cname
            for cname, cvals in df.items()
            if cvals.std() > 0
        ]

        view.plotly_box_multiple(
            mdata,
            **{
                "table": {
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
                        "tables": {
                            "value": [
                                f"{variable_name}.data"
                            ],
                            "sidebar": False
                        },
                        "filter_cols": {
                            "sidebar": False,
                            "type": {
                                "value": "index",
                                "sidebar": False
                            },
                            "tables": {
                                "value": [
                                    f"{variable_name}.data"
                                ],
                                "sidebar": False
                            },
                            "cname": {
                                "value": None,
                                "sidebar": False
                            },
                            "expr": {
                                "value": "in",
                                "sidebar": False
                            },
                            "value_enum": {
                                "value": values_to_show,
                                "sidebar": True
                            },
                            "value_str": {
                                "value": None,
                                "sidebar": False
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
                    "category": {
                        "enabled": {
                            "value": True,
                            "sidebar": False
                        },
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
                            "category": {
                                "sidebar": False,
                                "table": {
                                    "value": "Observation Metadata",
                                    "sidebar": False
                                },
                                "cname": {
                                    "value": "visitno",
                                    "sidebar": True
                                },
                                "label": {
                                    "value": "visitno",
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
                    }
                },
                "variable_options": {
                    "axis": {
                        "value": "Facet",
                        "sidebar": True
                    },
                    "log_values": {
                        "value": None,
                        "sidebar": True
                    },
                    "sort_by": {
                        "value": "Mean",
                        "sidebar": True
                    }
                },
                "category_options": {
                    "axis": {
                        "value": "Axis",
                        "sidebar": False
                    },
                    "sort_by": {
                        "value": "Mean",
                        "sidebar": False
                    }
                },
                "display_options": {
                    "ncols": {
                        "value": 4,
                        "sidebar": True
                    },
                    "outliers": {
                        "enabled": {
                            "value": True,
                            "sidebar": True
                        }
                    },
                    "title": {
                        "value": variable_name,
                        "sidebar": True
                    },
                    "var_label": {
                        "value": "Variable",
                        "sidebar": False
                    },
                    "val_label": {
                        "value": "Value",
                        "sidebar": False
                    },
                    "height": {
                        "value": 500,
                        "sidebar": False
                    },
                    "legend": {
                        "value": None,
                        "sidebar": False
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

    # Add the figures
    _add_figures(mdata)

    # Save the data
    _save_data(mdata)

    # Navigate to the figures
    st.switch_page("pages/view_all.py")


if __name__ == "__main__":
    main()
