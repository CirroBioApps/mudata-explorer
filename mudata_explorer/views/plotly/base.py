import pandas as pd
from mudata_explorer.base.view import View
import plotly.express as px


class Plotly(View):

    category = "Plotting"

    def fetch_dataframe(self, table_kw):
        """
        Helper method to fetch a dataframe from the params.
        This provides added functionality for parsing the color scales.
        """

        if self.params.get(f"{table_kw}.dataframe") is None:
            return None, {}

        data: pd.DataFrame = self.params[f"{table_kw}.dataframe"]
        colorscale = {}

        # Check each column to see if it is a color scale
        for col_kw in [
            kw[:-len(".scale")] for kw in self.params.keys()
            if kw.endswith(".scale") and kw.startswith(table_kw)
        ]:
            cname = col_kw[len(table_kw) + 1:]

            if self.params[f"{col_kw}.enabled"]:
                if self.params[f"{col_kw}.is_categorical"]:
                    data = data.assign(**{
                        cname: data[cname].apply(str)
                    })
                    colorscale = dict(
                        color_discrete_sequence=getattr(
                            px.colors.qualitative,
                            self.params[f"{col_kw}.scale"]
                        )
                    )
                else:
                    colorscale = dict(
                        color_continuous_scale=self.params[f"{col_kw}.scale"]
                    )
            else:
                colorscale = {}

        return data, colorscale
