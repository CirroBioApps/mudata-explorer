import pandas as pd
from mudata_explorer import app
from mudata_explorer.base.view import View
import seaborn as sns
from streamlit.delta_generator import DeltaGenerator


class Seaborn(View):

    category = "Plotting"


class Clustermap(Seaborn):

    type = "seaborn-clustermap"
    name = "Clustermap (Seaborn)"
    desc = "Display a heatmap with clustered rows and columnsusing Seaborn."
    schema = {
        "data": {
            "type": "dataframe",
            "select_columns": True,
            "query": "",
        },
        "z_score": {
            "type": "string",
            "label": "Z-Score Normalize",
            "help": "Normalize the data before clustering",
            "default": "None",
            "enum": ["None", "row", "column"]
        }
    }

    def display(self, container: DeltaGenerator):

        data: pd.DataFrame = self.params.get("data.dataframe")

        if data is None:
            return

        if data.shape[1] < 2:
            container.write("Please select at least two columns")
            return

        if self.params["z_score"] == "row":
            drop_rows = [
                obs for obs, vals in data.iterrows()
                if vals.std() == 0
            ]
            if len(drop_rows) > 0:
                if app.get_edit_views_flag():
                    container.write(
                        f"Dropping {len(drop_rows):,} rows with no variance."
                    )
                data = data.drop(drop_rows)

        elif self.params["z_score"] == "column":
            drop_cols = [
                col for col, vals in data.items()
                if vals.std() == 0
            ]
            if len(drop_cols) > 0:
                if app.get_edit_views_flag():
                    container.write(
                        "Dropping {:,} columns with no variance.".format(
                            len(drop_cols)
                        )
                    )
                data = data.drop(drop_cols, axis=1)

        g = sns.clustermap(
            data,
            cmap="viridis",
            figsize=(12, 12),
            z_score=(
                None if self.params["z_score"] == "None"
                else (
                    0 if self.params["z_score"] == "row"
                    else 1
                )
            )
        )
        container.pyplot(g)
