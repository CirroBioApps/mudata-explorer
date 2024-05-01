import anndata as ad
import json
import pandas as pd
import streamlit as st
from streamlit.delta_generator import DeltaGenerator
from mudata_explorer import app


def summarize_mdata(container: DeltaGenerator):

    container.write("**Current MuData**")

    mdata = app.get_mdata()
    if mdata is None:
        container.write("No data loaded.")
        return
    if len(mdata.mod) == 1 and "_blank" in mdata.mod:
        container.write("No data loaded.")
        return

    provenance = app.get_provenance(mdata)
    for mod_name, mod in mdata.mod.items():
        shape = mod.to_df().shape
        container.write(f" - {mod_name} ({shape[0]:,} observations x {shape[1]:,} features)") # noqa

        for slot in ["obsm", "obsp", "varm", "varp"]:

            for kw, df in getattr(mod, slot).items():
                container.write(f"   - {slot}[{kw}]: {df.shape[1]:,} cols")
                provenance_kw = app.format_provenance_key(mod_name, slot, kw)
                if provenance_kw in provenance:
                    event_str = json.dumps(provenance[provenance_kw])
                    container.write("   - " + event_str)


def display_table(container: DeltaGenerator):
    """Allow the user to display a table of the data."""

    mdata = app.get_mdata()
    if mdata is None or len(mdata.mod) == 1 and "_blank" in mdata.mod:
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
        df = pd.DataFrame(getattr(adata, attr)[kw], index=adata.obs.index)
    container.dataframe(df)


if __name__ == "__main__":
    app.setup_pages()

    st.write("#### Summarize")

    summarize_mdata(st.container())

    display_table(st.container())