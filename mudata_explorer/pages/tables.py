from mudata_explorer.base.slice import MuDataSlice
from mudata_explorer.app.sidebar import setup_sidebar
from mudata_explorer.helpers.join_kws import join_kws
from mudata_explorer.app.mdata import get_mdata_exists, setup_mdata, get_mdata, set_mdata, add_modality, list_modalities
from mudata_explorer.app.tables import get_table, build_mdata
from mudata_explorer.helpers.plotting import plot_mdata
import streamlit as st
from typing import Optional



def show_table(
    label: str,
    slot: str,
    mod: Optional[str] = None,
    attr: Optional[str] = None,
    editable=True,
    id="main"
):
    """Show the user a table and let them upload a new version."""

    kw = join_kws(mod, slot, attr)

    st.write(f"#### {label}")

    # Make a slice of the data
    slice = MuDataSlice(
        slot,
        modality=mod,
        axis=slot.startswith("var"),
        attr=attr
    )
    df = slice.dataframe(get_mdata(id=id, full=False))
    if isinstance(df, str):
        st.write(df)
        return

    if df is None or df.shape[0] == 0:
        if df is not None:
            st.write(df)
        st.write("No data currently available.")

    else:

        if df.shape[1] > 0:
            st.dataframe(df, hide_index=False)

        st.download_button(
            "Download (CSV)" if df.shape[1] > 0 else "Download Template (CSV)",
            df.to_csv(),
            f"{kw}.csv",
            mime="text/csv",
            help="Click here to download the data as a CSV file.",
            key=f"download_{kw}_csv"
        )

    # Skip the following if the table is not editable
    if not editable:
        return

    # Provide a button to add a new table
    if st.button(f"Add New {label} Table"):
        upload_csv_modal(slice, kw)


@st.dialog("Upload New Table", width="large")
def upload_csv_modal(slice: MuDataSlice, kw: str):

    # Let the user upload a new version
    new_df = get_table(
        st.container(),
        "Upload New Table (CSV, TSV, or XLSX)",
        keep_str=slice.slot == "obs",
        key=f"upload_{kw}"
    )

    if new_df is not None:

        st.write(
            "Read in {:,} rows x {:,} columns.".format(
                *new_df.shape
            )
        )

        if st.checkbox(
            f"Use All {new_df.shape[1]:,} Columns",
            value=True
        ):
            subset_columns = False
        else:
            subset_columns = st.multiselect(
                "Select Columns",
                new_df.columns.values,
                default=new_df.columns.values
            )

        if st.button("Add Table"):
            # Subset the columns if needed
            new_df = (
                new_df.reindex(columns=subset_columns)
                if subset_columns else new_df
            )

            # If there is no data, set it up
            if not get_mdata_exists():
                # A new dataset can only be set up with observation metadata
                # or a new set of observation data
                assert slice.slot == "obs", slice.slot
                setup_mdata()

            # Get the MuData object
            mdata = get_mdata(full=False)

            # Add the new data
            mdata = slice.write(
                mdata,
                new_df
            )

            set_mdata(mdata, full=False)
            st.rerun()


def show_modality(mod_name: str, id="main"):

    # Display the observation data, but don't let the user modify it
    show_table(
        label=f"{mod_name}: Observation Data",
        slot="X",
        mod=mod_name,
        editable=False,
        id=id
    )

    # Let the user modify the variable annotations
    show_table(
        label=f"{mod_name}: Variable Metadata",
        slot="var",
        mod=mod_name,
        id=id
    )

    # For any other tables which have been added
    for slot in ["varm", "varp", "obsm", "obsp"]:
        for attr in getattr(get_mdata(id=id, full=False).mod[mod_name], slot).keys():
            show_table(
                label=f"{mod_name}: {slot}[{attr}]",
                slot=slot,
                attr=attr,
                mod=mod_name
            )


@st.dialog("New Measurement Data", width="large")
def add_modality_modal(id="main"):

    # Let the user upload a new version
    new_df = get_table(
        st.container(),
        "Upload New Table",
        key="new_modality"
    )

    # Find the number of observations which overlap with the existing .obs
    if new_df is not None and get_mdata_exists():
        obs = get_mdata(id=id, full=False).obs
        overlap = set(obs.index).intersection(set(new_df.index))
        if len(overlap) > 0:
            st.write(
                f"There are {len(overlap):,} / {new_df.shape[0]:,} " +
                "observations which overlap with the existing data."
            )

    mod_name = st.text_input(
        "Measurement Name",
        help="Enter the name of the measurement.",
        value="Measurement"
    )
    if "." in mod_name:
        st.error("The measurement name cannot contain a period.")
        return

    if get_mdata_exists(id=id) and mod_name in list_modalities(id=id):
        st.write("A measurement with this name already exists.")
        return

    if new_df is not None and mod_name is not None and len(mod_name) > 0:

        st.write(
            "Read in {:,} rows x {:,} columns.".format(
                *new_df.shape
            )
        )

        if st.checkbox(
            f"Use All {new_df.shape[1]:,} Columns",
            value=True
        ):
            subset_columns = False
        else:
            subset_columns = st.multiselect(
                "Select Columns",
                new_df.columns.values,
                default=new_df.columns.values
            )

        if st.button(f"Add Table: {mod_name}"):

            if subset_columns:
                new_df = new_df.reindex(columns=subset_columns)

            # Add the new modality
            mdata = get_mdata(id=id, full=False)
            mdata = add_modality(
                mdata,
                mod_name,
                new_df
            )

            # Add to the global scope
            set_mdata(mdata, id=id, full=False)

            # Rerun the page
            st.rerun()


def run():
    """Display the page used to upload and modify data tables."""

    setup_sidebar()

    st.write("### MuData Contents")

    plot_mdata()

    # If no data exists
    if not get_mdata_exists():
        
        # Redirect to the load page
        st.switch_page("pages/load.py")

    # If data does exist
    else:

        # Show the metadata
        show_table(
            label="Observation Metadata",
            slot="obs"
        )

        # Show each of the measurement modalities
        for mod_name in list_modalities():
            show_modality(mod_name)
            st.write("---")
