from mudata_explorer.base.view import View
from streamlit.delta_generator import DeltaGenerator


class SummarizeMuData(View):

    type = "summarize-mudata"
    name = "Summarize MuData"
    help_text = "Print a short description of the data available."
    category = "Summary"
    defaults = {}

    def display(self, container: DeltaGenerator):

        mdata = self.get_mdata()

        if (
            mdata is None or
            mdata.shape[0] == 0 or
            (mdata.shape[0] == 1 and len(mdata.mod) == 1 and "_blank" in mdata.mod)
        ):
            container.write("No MuData object available.")
            return

        container.write(f"{mdata.shape[0]:,} observations")
        container.write(f"{mdata.obs.shape[1]:,} metadata annotations")
        for key, adata in mdata.mod.items():
            container.write(f" - {key}: {adata.shape[1]:,} measurements")
            for kw in adata.obsm.keys():
                container.write(f"   - obsm[{kw}]: {adata.obsm[kw].shape[1]:,} observations")
            for kw in adata.varm.keys():
                container.write(f"   - varm[{kw}]: {adata.varm[kw].shape[1]:,} observations")
