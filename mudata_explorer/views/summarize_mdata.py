from mudata_explorer.base.view import View
from streamlit.delta_generator import DeltaGenerator


class SummarizeMuData(View):

    type = "summarize-mudata"
    name = "Summarize MuData"
    desc = "Print a short description of the data available."
    categories = ["Summary"]
    defaults = {}

    def display(self, container: DeltaGenerator):

        mdata = self.get_mdata()

        if mdata is None or mdata.shape[0] == 0:
            container.write("No MuData object available.")
            return

        container.write(f"{mdata.shape[0]:,} observations")
        container.write(f"{mdata.obs.shape[1]:,} metadata annotations")
        for key, adata in mdata.mod.items():
            container.write(f" - {key}: {adata.shape[1]:,} measurements")
