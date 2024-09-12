import anndata as ad
import pandas as pd
import streamlit as st
from streamlit.delta_generator import DeltaGenerator
from mudata_explorer.helpers import join_kws
from mudata_explorer.app.mdata import get_mdata, get_mdata_exists, get_provenance
from mudata_explorer.app.sidebar import setup_sidebar


def summarize_mdata(container: DeltaGenerator, id="main"):

    container.write("**Current MuData**")

    if not get_mdata_exists(id=id):
        container.write("No data loaded.")
        return

    mdata = get_mdata(id=id, full=False)

    provenance = get_provenance(id=id)

    for mod_name, mod in mdata.mod.items():
        shape = mod.to_df().shape
        container.write(f" - {mod_name} ({shape[0]:,} observations x {shape[1]:,} features)") # noqa
        show_provenance(mod_name, "X", None, provenance, container)

        for slot in ["obsm", "obsp", "varm", "varp"]:

            for kw, df in getattr(mod, slot).items():
                container.write(f"   - {slot}[{kw}]: {df.shape[1]:,} cols")
                show_provenance(mod_name, slot, kw, provenance, container)


def show_provenance(
    mod_name: str,
    slot: str,
    kw: str,
    provenance: dict,
    container: DeltaGenerator
):
    provenance_key = join_kws(mod_name, slot, kw)
    if provenance.get(provenance_key) is not None:
        with container.expander("Show provenance"):
            st.write(provenance.get(provenance_key))


def display_table(container: DeltaGenerator, id="main"):
    """Allow the user to display a table of the data."""

    if not get_mdata_exists(id=id):
        return

    mdata = get_mdata(id=id, full=False)

    container.write("#### Display Table")

    # Select the modality to display
    modality = container.selectbox(
        "Select modality",
        list(mdata.mod.keys())
    )

    # Get the data for the selected modality
    adata: ad.AnnData = mdata.mod[modality]

    # Let the user pick which table to show
    table_options = ["metadata", "data"]
    for attr in ["obsm", "obsp", "varm", "varp"]:
        for kw in getattr(adata, attr).keys():
            table_options.append(f"{attr}[{kw}]")
    table = container.selectbox(
        "Select table",
        table_options
    )

    if table == "metadata":
        df = adata.obs
    elif table == "data":
        df = adata.to_df()
    else:
        attr, kw = table.split("[")
        kw = kw[:-1]
        df = pd.DataFrame(getattr(adata, attr)[kw], index=adata.obs_names)
    container.dataframe(df)

    container.download_button(
        "Download table as CSV",
        df.to_csv(),
        f"{modality}_{table}.csv",
        "text/csv",  # Add the MIME type for CSV
    )


def run():
    setup_sidebar("summarize")

    st.write("#### Summarize")

    summarize_mdata(st.container())

    display_table(st.container())
