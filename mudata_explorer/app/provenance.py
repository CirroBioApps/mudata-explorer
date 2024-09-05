import streamlit as st
from mudata_explorer.base.slice import MuDataSlice
from mudata_explorer.app.mdata import query_provenance
from mudata_explorer.app.process import process_sdk_snippet
from streamlit.delta_generator import DeltaGenerator
from plotly import io


def show_provenance(loc: MuDataSlice, container: DeltaGenerator, id="main"):

    prov = query_provenance(loc, id=id)
    if prov is not None:
        with container.expander(
            f"**Provenance: '{loc.address}'**"
        ):
            st.write(f"**Analysis:** {prov['process']}")
            st.write(f"**Timestamp:** {prov['timestamp']}")
            st.write("**SDK Configuration:**")
            st.code(process_sdk_snippet(prov))

            if isinstance(prov.get("figures"), list):
                if len(prov.get("figures")) > 0:
                    st.write("Supporting Figures")
                for fig_json in prov["figures"]:
                    fig = io.from_json(fig_json)
                    st.plotly_chart(fig)
        return True
    return False

