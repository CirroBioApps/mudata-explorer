from mudata_explorer import app
from mudata_explorer.helpers.join_kws import join_kws
from typing import Union
import pandas as pd
import plotly.express as px
from streamlit.delta_generator import DeltaGenerator
import streamlit as st


class MuDataAppHelpers:
    params: dict = {}
    schema: dict
    orientation_schema: dict = {
        "orientation": {
            "type": "string",
            "enum": ["observations", "variables"],
            "default": "observations",
            "label": "Orientation - Assign Outputs To:",
            "help": "Select whether to use observations or variables."
        }
    }

    def param(self, *kws, default=None):
        return self.params.get(join_kws(*kws), default)

    def get_schema_defaults(self, schema: dict, prefix=None):
        """Yield the default values for each item in the schema."""

        # Orientation is an essential parameter
        yield "orientation", self.orientation_schema["orientation"]["default"]

        for key, elem in schema.items():
            assert isinstance(elem, dict), f"Expected dict, got {type(elem)}"
            assert "type" in elem, f"Missing 'type' key in schema: {elem}"

            if elem["type"] == "object":
                yield from self.get_schema_defaults(
                    elem["properties"],
                    join_kws(prefix, key)
                )

            elif elem["type"] in ["string", "number", "float", "boolean"]:
                yield (
                    join_kws(prefix, key),
                    elem.get("default", None)
                )

            elif elem["type"] == "dataframe":
                # If column selection is not enabled,
                # the user will select >= 1 tables
                if len(elem.get("columns", {})) == 0:
                    yield (
                        join_kws(prefix, key, "tables"),
                        elem.get("tables", [])
                    )
                # If the user can filter down the columns of interest
                if elem.get("select_columns", False):
                    yield (
                        join_kws(prefix, key, "selected_columns"),
                        []
                    )
                # Add any columns specified in the schema
                for col_kw, col_elem in elem.get("columns", {}).items():
                    for kw in ["modality", "table", "cname", "label"]:
                        yield (
                            join_kws(prefix, key, col_kw, kw),
                            col_elem.get(kw, None)
                        )
                    if col_elem.get("optional", False):
                        yield (
                            join_kws(prefix, key, col_kw, "enabled"),
                            True
                        )
                    if col_elem.get("continuous_scale", False):
                        yield (
                            join_kws(prefix, key, col_kw, "continuous_scale"),
                            "Viridis"
                        )
                    if col_elem.get("discrete_sequence", False):
                        yield (
                            join_kws(prefix, key, col_kw, "discrete_sequence"),
                            "Plotly"
                        )
                # If the user is allowed to filter the data
                if elem.get("query", True):
                    for attr in [
                        "modality",
                        "table",
                        "cname",
                        "expr",
                        "value"
                    ]:
                        yield (
                            join_kws(prefix, key, "query", attr),
                            ""
                        )

    def input_value_kwargs(self, kw, copy_to=None):
        """Each input value element will be populated with default kwargs."""
        return dict(
            value=self.params[kw],
            key=self.param_key(kw),
            on_change=self.input_value_change,
            args=(kw,),
            kwargs=dict(copy_to=copy_to)
        )

    def update_view_param(self, kw, value):
        # Get the MuData object
        mdata = app.get_mdata()

        # Modify the value of this param for this view
        mdata.uns["mudata-explorer-views"][self.ix]["params"][kw] = value

        # Save the MuData object
        app.set_mdata(mdata)

        # Also update the params object
        self.params[kw] = value

    def input_selectbox_kwargs(self, kw, options: list, copy_to=None):
        """Populate the selectbox element with default kwargs."""
        if self.params[kw] not in options:
            index = 0
        else:
            index = options.index(self.params[kw])

        if self.params[kw] != options[index]:
            self.update_view_param(kw, options[index])

        return dict(
            index=index,
            key=self.param_key(kw),
            on_change=self.input_value_change,
            args=(kw,),
            kwargs=dict(copy_to=copy_to)
        )

    def input_multiselect_kwargs(self, kw, options, copy_to=None):
        """Populate the multiselect element with default kwargs."""

        default = [
            val for val in self.params[kw]
            if val in options
        ]

        return dict(
            key=self.param_key(kw),
            default=default,
            on_change=self.input_value_change,
            args=(kw,),
            kwargs=dict(copy_to=copy_to)
        )

    def input_value_change(self, kw, copy_to=None):
        # Get the value provided by the user
        value = st.session_state[self.param_key(kw)]

        # If this is different than the params
        if value != self.params[kw]:

            # Update the view in the mdata object
            self.update_view_param(kw, value)

            # If a copy_to kwarg was provided
            if copy_to is not None:
                # Modify that value as well
                self.update_view_param(copy_to, value)

    def get_data(self, container: DeltaGenerator):

        container.write("##### Inputs")

        # 'orientation' is a special-case parameter that is used to
        # determine whether the user is selecting observations or variables
        msg = "The 'orientation' parameter is reserved."
        assert "orientation" not in self.schema.keys(), msg
        self.render_form(
            container,
            self.orientation_schema
        )

        # Parse the form schema of the object
        return self.render_form(container, self.schema)

    def render_form(
        self,
        container: DeltaGenerator,
        schema: dict,
        prefix: str = None
    ):
        # Get the global settings
        settings = app.get_settings()

        # Iterate over the form defined for this view
        for key, elem in schema.items():

            prefix_key = join_kws(prefix, key)

            if elem["type"] == "dataframe":

                self.render_dataframe(prefix, key, elem, container)

            elif elem["type"] == "object":

                if "label" in elem and settings["editable"]:
                    container.write(f"#### {elem['label']}")

                self.render_form(
                    container,
                    elem["properties"],
                    prefix_key
                )

            elif elem["type"] == "string":

                if settings["editable"]:
                    if elem.get("enum") is not None:
                        self.params[prefix_key] = container.selectbox(
                            key if elem.get("label") is None else elem["label"],
                            elem["enum"],
                            help=elem.get("help"),
                            **self.input_selectbox_kwargs(
                                prefix_key,
                                elem["enum"]
                            )
                        )

                    elif elem.get("multiline", False):
                        self.params[prefix_key] = container.text_area(
                            (
                                key
                                if elem.get("label") is None
                                else elem["label"]
                            ),
                            help=elem.get("help"),
                            **self.input_value_kwargs(prefix_key)
                        )

                    else:
                        self.params[prefix_key] = container.text_input(
                            (
                                key
                                if elem.get("label") is None
                                else elem["label"]
                            ),
                            help=elem.get("help"),
                            **self.input_value_kwargs(prefix_key)
                        )

            elif elem["type"] in ["number", "float"]:

                if settings["editable"]:
                    self.params[prefix_key] = container.number_input(
                        key if elem.get("label") is None else elem["label"],
                        help=elem.get("help"),
                        **self.input_value_kwargs(prefix_key)
                    )

            elif elem["type"] == "boolean":

                if settings["editable"]:
                    self.params[prefix_key] = container.checkbox(
                        key if elem.get("label") is None else elem["label"],
                        help=elem.get("help"),
                        **self.input_value_kwargs(prefix_key)
                    )

            else:
                raise Exception(f"Unsupported type: {elem['type']}")

    def render_dataframe(
        self,
        prefix: Union[None, str],
        key: str,
        elem: dict,
        container: DeltaGenerator
    ):
        # Get the global settings
        settings = app.get_settings()

        mdata = app.get_mdata()
        if mdata is None or mdata.shape[0] == 0:
            container.write("No MuData object available.")
            return

        # The user can either:
        # 1. Select columns to use
        # 2. Select 1 or more tables in their entirety
        # 3. Select specific column from 1 or more tables

        # If 'columns' were specified
        if len(elem.get("columns", {})) > 0:

            # Iterate over the columns in the schema and prompt the user
            # for which data to supply for each
            for col_kw, col_elem in elem.get("columns", {}).items():

                self.render_dataframe_column(
                    prefix,
                    key,
                    col_kw,
                    col_elem,
                    container
                )

            df = self.build_dataframe(join_kws(prefix, key), elem["columns"])

        # If 'columns' was not specified
        else:

            tables_kw = join_kws(prefix, key, "tables")

            # Let the user select one or more tables
            if settings["editable"]:

                all_tables = app.tree_tables(self.params["orientation"])
                container.multiselect(
                    "Select table(s)",
                    all_tables,
                    **self.input_multiselect_kwargs(
                        tables_kw,
                        all_tables
                    )
                )

            # Make a DataFrame with the selected table(s)
            selected_tables = self.params.get(tables_kw, [])
            if len(selected_tables) == 0:
                container.write("No tables selected.")
                return
            else:
                df = pd.concat(
                    [
                        app.get_dataframe_table(*table_path.split(".", 1))
                        for table_path in self.params.get(tables_kw, [])
                    ],
                    axis=1
                )

            # If the orientation is set to variables, transpose the DataFrame
            if self.params["orientation"] == "variables":
                df = df.T

        # Let the user optionally filter rows
        if elem.get("query", True):
            filtered_obs = self.render_query(prefix, key, container)
            if filtered_obs is not None:
                df = df.loc[filtered_obs]
                if settings["editable"]:
                    container.write(f"Filtered to {df.shape[0]:,} samples.")

        # If the user has the option to select specific columns
        if elem.get("select_columns", False):
            # Get the list of available columns
            avail_columns = list(df.columns.values)

            # Get the columns selected by the user
            selected_columns_kw = join_kws(prefix, key, "selected_columns")
            selected_columns = self.param(selected_columns_kw, default=[])

            # Make sure that the selected columns are valid
            if any([cname not in avail_columns for cname in selected_columns]):
                selected_columns = [
                    cname
                    for cname in selected_columns
                    if cname in avail_columns
                ]
                self.update_view_param(selected_columns_kw, selected_columns)

            # If no columns are selected
            if len(selected_columns) == 0:
                # Select all of the coluns
                self.update_view_param(selected_columns_kw, avail_columns)

            # Let the user select the columns to use
            if settings["editable"]:
                # Check to see if the user wants to select all of the columns
                if container.checkbox(
                    "Use all columns",
                    set(selected_columns) == set(avail_columns)
                ):
                    self.update_view_param(selected_columns_kw, avail_columns)

                # Otherwise, present a list of columns to select from
                else:
                    container.multiselect(
                        "Select columns",
                        avail_columns,
                        **self.input_multiselect_kwargs(
                            selected_columns_kw,
                            avail_columns
                        )
                    )

            # Filter the DataFrame to the selected columns
            df = df[selected_columns]

        if settings["editable"] and df.shape[0] == 0:
            container.write("No data available.")
            return

        self.params[join_kws(prefix, key, "dataframe")] = df

    def render_dataframe_column(
        self,
        prefix: Union[None, str],
        key: str,
        col_kw: str,
        col_elem: dict,
        container: DeltaGenerator
    ):

        # Get the global settings
        settings = app.get_settings()

        # Set the default values

        # Modality selection
        all_modalities = app.list_modalities()
        mod_kw = join_kws(prefix, key, col_kw, "modality")
        if self.params.get(mod_kw) is None:
            self.update_view_param(
                mod_kw,
                all_modalities[0]
            )

        # Table selection
        table_kw = join_kws(prefix, key, col_kw, "table")
        if self.params.get(table_kw) is None:
            self.update_view_param(
                table_kw,
                "data"
            )

        # Column label selection
        label_kw = join_kws(prefix, key, col_kw, "label")
        if self.params.get(label_kw) is None:
            self.update_view_param(
                label_kw,
                "Label"
            )

        # Column name selection
        cname_kw = join_kws(prefix, key, col_kw, "cname")
        if self.params.get(cname_kw) is None:
            cname_val = app.list_cnames(
                self.params[mod_kw],
                self.params[table_kw],
                orientation=self.params["orientation"]
            )[0]

            self.update_view_param(cname_kw, cname_val)
            self.update_view_param(label_kw, cname_val)

        # If the views are editable
        if settings["editable"]:

            # Print the column name
            container.write(f"#### {col_elem.get('label', col_kw)}")

            # Optional column selection
            enabled_kw = join_kws(prefix, key, col_kw, "enabled")
            if col_elem.get("optional", False):
                container.checkbox(
                    "Enabled",
                    **self.input_value_kwargs(enabled_kw)
                )

            # If the column is enabled
            if self.params.get(enabled_kw, True):

                # Make three columns for the modality, table, and column name
                cols = container.columns([1, 1, 1, 1])

                # Select the modality
                cols[0].selectbox(
                    "Modality",
                    all_modalities,
                    **self.input_selectbox_kwargs(mod_kw, all_modalities)
                )

                # Get the list of tables available for this modality
                all_tables = app.list_tables(
                    self.params[mod_kw],
                    self.params["orientation"]
                )

                # Select the table of interest
                cols[1].selectbox(
                    "Table",
                    all_tables,
                    **self.input_selectbox_kwargs(table_kw, all_tables)
                )

                # Get the list of possible columns
                all_cnames = app.list_cnames(
                    self.params[mod_kw],
                    self.params[table_kw]
                )

                # Select the column name
                cols[2].selectbox(
                    "Column",
                    all_cnames,
                    **self.input_selectbox_kwargs(
                        cname_kw,
                        all_cnames,
                        copy_to=label_kw
                    )
                )

                # Input the column label
                cols[3].text_input(
                    "Label",
                    **self.input_value_kwargs(label_kw)
                )

                # Color options
                if col_elem.get("continuous_scale", False):

                    container.selectbox(
                        "Select color scale",
                        px.colors.named_colorscales(),
                        **self.input_selectbox_kwargs(
                            join_kws(prefix, key, col_kw, "continuous_scale"),
                            px.colors.named_colorscales()
                        )
                    )

                if col_elem.get("discrete_sequence", False):

                    container.selectbox(
                        "Select color sequence",
                        dir(px.colors.qualitative),
                        **self.input_selectbox_kwargs(
                            join_kws(prefix, key, col_kw, "discrete_sequence"),
                            dir(px.colors.qualitative)
                        )
                    )

    def build_dataframe(self, key: str, columns: dict):
        # Get the information for each column
        return pd.DataFrame({
            col_kw: app.get_dataframe_column(
                **{
                    kw: self.param(key, col_kw, kw)
                    for kw in ["modality", "table", "cname"]
                }
            )
            for col_kw, col_elem in columns.items()
            if (
                col_elem.get("optional", False) is False
                or
                self.param(key, col_kw, "enabled")
            )
        })

    def render_query(self, prefix: str, key: str, container: DeltaGenerator):

        # Get the global settings
        settings = app.get_settings()

        # Set the default values

        # Modality selection
        all_modalities = app.list_modalities()
        mod_kw = join_kws(prefix, key, "query", "modality")
        if self.params.get(mod_kw) is None:
            self.update_view_param(
                mod_kw,
                all_modalities[0]
            )

        # Table selection
        table_kw = join_kws(prefix, key, "query", "table")
        if self.params.get(table_kw) is None:
            self.update_view_param(
                table_kw,
                "data"
            )

        # Column name selection
        cname_kw = join_kws(prefix, key, "query", "cname")
        if self.params.get(cname_kw) is None:

            self.update_view_param(
                cname_kw,
                app.list_cnames(
                    self.params[mod_kw],
                    self.params[table_kw]
                )[0]
            )

        # Boolean operator selection
        expr_kw = join_kws(prefix, key, "query", "expr")
        if self.params.get(expr_kw) is None:

            self.update_view_param(
                expr_kw,
                ">="
            )

        # Boolean value selection
        value_kw = join_kws(prefix, key, "query", "value")
        if self.params.get(value_kw) is None:

            self.update_view_param(value_kw, "")

        # If the views are editable
        if settings["editable"]:

            # Print the column name
            container.write("#### Filter samples")

            # Make three columns for the modality, table, and column name
            cols = container.columns([1, 1, 1, 1, 1])

            # Select the modality
            cols[0].selectbox(
                "Modality",
                all_modalities,
                **self.input_selectbox_kwargs(mod_kw, all_modalities)
            )

            # Get the list of tables available for this modality
            all_tables = app.list_tables(
                self.params[mod_kw],
                self.params["orientation"]
            )

            # Select the table of interest
            cols[1].selectbox(
                "Table",
                all_tables,
                **self.input_selectbox_kwargs(table_kw, all_tables)
            )

            # Get the list of possible columns
            all_cnames = app.list_cnames(
                self.params[mod_kw],
                self.params[table_kw]
            )

            # Select the column name
            cols[2].selectbox(
                "Column",
                all_cnames,
                **self.input_selectbox_kwargs(
                    cname_kw,
                    all_cnames
                )
            )

            # Input the boolean operator
            cols[3].selectbox(
                "Operator",
                [">=", "<=", "==", "!=", ">", "<"],
                **self.input_selectbox_kwargs(expr_kw, [">=", "<=", "==", "!=", ">", "<"])
            )

            # Input the boolean value
            cols[4].text_input(
                "Value",
                help="Enter a value to filter samples on.",
                **self.input_value_kwargs(value_kw)
            )

        # Get the values for the query
        query = {
            "modality": self.param(key, "query", "modality"),
            "table": self.param(key, "query", "table"),
            "cname": self.param(key, "query", "cname"),
            "expr": self.param(key, "query", "expr"),
            "value": self.param(key, "query", "value")
        }

        # If no value is provided
        if query['value'] is None or len(query['value']) == 0:
            if settings["editable"]:
                container.write("Provide a value to filter samples.")
            return

        # Get the table
        table = app.get_dataframe_table(
            query["modality"],
            query["table"]
        )
        if table is None:
            if settings["editable"]:
                container.write("No data available for filtering.")
            return

        # Apply the filter
        try:
            table = table.query(f"{query['cname']} {query['expr']} {query['value']}")
        except Exception as e:
            container.write("Error while filtering")
            container.exception(e)
            return

        # If no values are returned
        if table.shape[0] == 0:
            container.write("No samples match the filter criteria.")
            return

        return table.index
