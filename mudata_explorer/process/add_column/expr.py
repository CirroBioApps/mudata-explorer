import pandas as pd
from mudata_explorer.base.process import Process


class RunAddColumnExpression(Process):

    type = "add_column_expr"
    name = "Add Column: Expression"
    help_text = """
Add a column which computes an arbitrary expression on the input data.

For example, if you have a data table with values like:

| length | width | height |
|--------|-------|--------|
| 10     | 20    | 30     |
| 15     | 25    | 35     |
| 20     | 30    | 40     |

Then you could add a column which computes the volume:

`length * width * height`

Which would result in:

| length | width | height | volume |
|--------|-------|--------|--------|
| 10     | 20    | 30     | 6000   |
| 15     | 25    | 35     | 13125  |
| 20     | 30    | 40     | 24000  |

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
        "expression": {
            "type": "string",
            "label": "Expression",
            "help": "Enter the expression to apply to the data",
            "placeholder": "length * width * height",
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
            "desc": "Results of the user-supplied expression.",
            "modality": "table.data.tables",
            "axis": "table.data.axis",
            "attr": "outputs.dest_key"
        }
    }

    def execute(self):

        df: pd.DataFrame = self.params["table.data.dataframe"]

        # Apply the formula to the table
        res = df.eval(self.params["expression"])
        assert isinstance(res, pd.Series), df.eval(self.params["formula"])
        
        self.save_results("results", res)
