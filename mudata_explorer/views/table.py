import pandas as pd
from mudata_explorer.base.view import View
import streamlit as st


class Table(View):

    category = "Summary"
    type = "table"
    name = "Table"
    help_text = "Show a table of data."
    schema = {
        "data": {
            "type": "object",
            "label": "Data Table",
            "properties": {
                "table": {
                    "type": "dataframe",
                    "label": "Data",
                    "select_columns": True,
                    "query": "",
                    "sidebar": True
                }
            }
        },
        "options": {
            "type": "object",
            "label": "Options",
            "properties": {
                "sort": {
                    "type": "dataframe",
                    "columns": {
                        "sort_by": {"label": "sort_by"}
                    },
                    "query": True,
                    "sidebar": True
                }
            }
        }
    }

    def display(self):

        data: pd.DataFrame = self.params.get("data.table.dataframe")

        if data is None or data.shape[1] < 1:
            st.write("Please select at least one column")
            return

        # Sort by the column specified in the sort_by column
        sort_by: pd.DataFrame = self.params.get("options.sort.dataframe")
        if sort_by is None or sort_by.shape[1] < 1:
            st.write("Please select a column to sort by")
            return

        # If the sort column is transposed, fix it
        if len(sort_by.columns.intersection(data.index)) > len(sort_by.index.intersection(data.index)):
            sort_by = sort_by.T

        # Sort the data by the column specified in the sort_by column
        data = data.reindex(sort_by.sort_values(by="sort_by").index)

        st.dataframe(data)
