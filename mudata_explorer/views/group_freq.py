import pandas as pd
from mudata_explorer.base.view import View
import streamlit as st


class GroupFreq(View):

    category = "Summary"
    type = "group-freq"
    name = "Group Frequency"
    help_text = """Show a table of the frequency of every combination
    of unique values across any number of selected columns.

This display is most effective for a small number of categorical values
in which the frequency of each unique combination is of interest.

Example input:

| Column A | Column B | Column C |
|----------|----------|----------|
| A        | X        | Y        |
| A        | Y        | Y        |
| B        | X        | Y        |
| B        | X        | Y        |
| B        | Y        | Y        |
| B        | Y        | Y        |

Example output:

| Column A | Column B | Column C | Frequency |
|----------|----------|----------|-----------|
| A        | X        | Y        | 1         |
| A        | Y        | Y        | 1         |
| B        | X        | Y        | 2         |
| B        | Y        | Y        | 2         |

    """
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
                    "type": "string",
                    "label": "Sort By",
                    "default": "Frequency",
                    "enum": ["Frequency", "Groupings"],
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

        freq = (
            data
            .groupby(list(data.columns))
            .size()
            .reset_index(name="Frequency")
            .query("Frequency > 0")
        )

        if self.params["options.sort"] == "Frequency":
            freq = freq.sort_values("Frequency", ascending=False)
        else:
            freq = freq.sort_values(list(data.columns))

        st.dataframe(freq, hide_index=True)
