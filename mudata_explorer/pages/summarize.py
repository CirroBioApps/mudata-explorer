import anndata as ad
import pandas as pd
import streamlit as st
from streamlit.delta_generator import DeltaGenerator
from mudata_explorer import app
from mudata_explorer.helpers import join_kws


def summarize_mdata(container: DeltaGenerator):

    container.write("**Current MuData**")

    mdata = app.get_mdata()
    if not app.has_mdata():
        container.write("No data loaded.")
        return

    provenance = app.get_provenance()

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


def display_table(container: DeltaGenerator):
    """Allow the user to display a table of the data."""

    mdata = app.get_mdata()
    if not app.has_mdata():
        return

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
    app.setup_sidebar()

    st.write("#### Summarize")

    summarize_mdata(st.container())

    display_table(st.container())
