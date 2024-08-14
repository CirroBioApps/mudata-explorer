from mudata_explorer.base.view import View
from mudata_explorer.app.mdata import query_provenance
from streamlit.delta_generator import DeltaGenerator
from mudata_explorer.base.slice import MuDataSlice
from plotly import io as pio


class SupportingFigure(View):

    category = "Supporting Figure"
    type = "supporting_figure"
    name = "Supporting Figure"
    help_text = "Display a figure which was generated during an analysis process."
    schema = {
        "input": {
            "type": "supporting_figure"
        }
    }

    def display(self, container: DeltaGenerator):

        # Get the location and index position of the figure
        loc_ix: str = self.params.get("input")

        if loc_ix is None:
            return

        # Split the location and the index
        loc, ix = loc_ix.rsplit(":", 1)
        ix = int(ix)

        # Get the provenance from the specified location
        provenance = query_provenance(
            MuDataSlice.hydrate(loc)
        )

        if provenance is None:
            container.write("Data not found.")

        # Get the figure from the provenance
        figures = provenance.get("figures", [])

        if len(figures) < (ix + 1):
            container.write("Figure not found.")

        # Show the figure
        container.plotly_chart(
            pio.from_json(figures[ix])
        )
