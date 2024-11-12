from typing import Optional
import pandas as pd
from mudata_explorer.base.view import View
import streamlit as st
from datastory.datastory import DataStory


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
        },
        "display_options": {
            "type": "object",
            "label": "Display Options",
            "properties": {
                "title": {
                    "type": "string",
                    "label": "Title",
                    "default": "",
                    "sidebar": True
                },
                "legend": {
                    "type": "string",
                    "label": "Legend",
                    "multiline": True
                }
            }
        }
    }

    def _format_table(self) -> Optional[pd.DataFrame]:

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

        return data

    def display(self):

        data = self._format_table()
        if data is None:
            return

        title = self.params["display_options.title"]
        if title:
            st.markdown(f"**{title}**")

        st.dataframe(data)

        legend = self.params["display_options.legend"]
        if legend:
            st.markdown(legend)

    def to_datastory(self, ds: DataStory):
        """
        Convert the view to a DataStory object.
        """
        
        self.params = self.form.dump()

        data = self._format_table()
        if data is None:
            return

        title = self.params["display_options.title"]
        if title:
            ds.add_markdown(
                f"**{title}**",
                style={"flex-basis": "80%"},
                section_style={"justify-content": "center", "align-items": "center"}
            )

        ds.add_dataframe(
            data,
            style={"flex-basis": "80%"},
            section_style={"justify-content": "center", "align-items": "center"}
        )
