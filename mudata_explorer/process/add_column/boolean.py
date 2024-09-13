import pandas as pd
from mudata_explorer.base.process import Process


class RunAddColumnBoolean(Process):

    type = "add_column_boolean"
    name = "Add Column: Boolean"
    help_text = """
Add a column which contains a boolean (true/false) value based on a condition.

For example, if you have metadata on a set of observations which includes:

| age | city |
|-----| ---- |
| 25  | NYC  |
| 30  | NYC  |
| 40  | LA   |
| 50  | LA   |

Then you could add a column which indicates if the age is greater than 25:

`age > 25`

Which would result in:

| age | city | age over 25 |
|-----| ---- | ----------- |
| 25  | NYC  | False       |
| 30  | NYC  | True        |
| 40  | LA   | True        |
| 50  | LA   | True        |

(note that the column name can be customized)

If instead you wanted to add a column which indicates if the age is over 25 and the city is NYC:

`age > 25 and city == 'NYC'`

Which would result in:

| age | city | age over 25 and city is NYC |
|-----| ---- | --------------------------- |
| 25  | NYC  | False                       |
| 30  | NYC  | True                        |
| 40  | LA   | False                       |
| 50  | LA   | False                       |

    """ # noqa

    category = "Add Column"
    schema = {
        "table": {
            "type": "object",
            "label": "Data Table",
            "properties": {
                "data": {
                    "type": "dataframe",
                    "label": "Data",
                    "help": "Select the data which will be used apply the formula", # noqa
                    "select_columns": True,
                    "query": "",
                }
            }
        },
        "formula": {
            "type": "string",
            "label": "Formula",
            "help": "Enter the formula to apply to the data",
            "placeholder": "age > 25",
        },
        "outputs": {
            "type": "object",
            "label": "Outputs",
            "properties": {
                "dest_key": {
                    "type": "string",
                    "default": "new_column",
                    "label": "Label to use for results",
                    "help": """
                    Key to use when saving the output
                    """
                }
            }
        }
    }
    outputs = {
        "results": {
            "type": pd.Series,
            "label": "Boolean Column",
            "desc": "Results of the user-supplied formula.",
            "modality": "table.data.tables",
            "axis": "table.data.axis",
            "attr": "outputs.dest_key"
        }
    }

    def execute(self):

        df: pd.DataFrame = self.params["table.data.dataframe"]

        # Apply the formula to the table
        res = df.eval(self.params["formula"])
        assert isinstance(res, pd.Series), df.eval(self.params["formula"])
        
        self.save_results("results", res)
