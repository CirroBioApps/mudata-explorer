import json
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


if __name__ == "__main__":
    app.setup_pages()

    st.write("#### Settings")

    summarize_mdata(st.container())
