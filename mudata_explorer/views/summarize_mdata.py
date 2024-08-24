from mudata_explorer.base.view import View
import streamlit as st


class SummarizeMuData(View):

    type = "summarize-mudata"
    name = "Summarize MuData"
    help_text = "Print a short description of the data available."
    category = "Summary"
    defaults = {}

    def display(self):

        mdata = self.get_mdata()

        if (
            mdata is None or
            mdata.shape[0] == 0 or
            (mdata.shape[0] == 1 and len(mdata.mod) == 1 and "_blank" in mdata.mod)
        ):
            st.write("No MuData object available.")
            return

        st.write(f"{mdata.shape[0]:,} observations")
        st.write(f"{mdata.obs.shape[1]:,} metadata annotations")
        for key, adata in mdata.mod.items():
            st.write(f" - {key}: {adata.shape[1]:,} measurements")
            for kw in adata.obsm.keys():
                st.write(f"   - obsm[{kw}]: {adata.obsm[kw].shape[1]:,} observations")
            for kw in adata.varm.keys():
                st.write(f"   - varm[{kw}]: {adata.varm[kw].shape[1]:,} observations")
