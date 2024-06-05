from mudata_explorer import app
from mudata_explorer.helpers.join_kws import join_kws
from typing import Dict, List, Optional, Union
import pandas as pd
import plotly.express as px
from streamlit.delta_generator import DeltaGenerator
import streamlit as st


class MuDataAppHelpers:

    ix: int
    params: dict = {}
    schema: dict
    type: str
    name: str
    desc: str
    category: str
    # Flag used to indicate whether all parameters have been provided
    params_complete: bool
    # Flag used to indicate whether the parameters are editable
    params_editable: bool

    def param(self, *kws, default=None):
        return self.params.get(join_kws(*kws), default)

    def get_schema_defaults(self, schema: dict, prefix=None):
        """Yield the default values for each item in the schema."""

        for key, elem in schema.items():
            assert isinstance(elem, dict), f"Expected dict, got {type(elem)}"
            assert "type" in elem, f"Missing 'type' key in schema: {elem}"

            if elem["type"] == "object":
                yield from self.get_schema_defaults(
                    elem["properties"],
                    join_kws(prefix, key)
                )

            elif elem["type"] in [
                "string",
                "float",
                "boolean",
                "integer",
                "supporting_figure"
            ]:
                yield (
                    join_kws(prefix, key),
                    elem.get("default", None)
                )
                if elem.get("optional", False):
                    yield (
                        join_kws(prefix, key, "enabled"),
                        True
                    )

            elif elem["type"] == "dataframe":
                # Each dataframe may be oriented to the obs or var
                yield (
                    join_kws(prefix, key, "axis"),
                    elem.get("axis", 0)
                )
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
                    for kw in ["table", "cname", "label"]:
                        yield (
                            join_kws(prefix, key, col_kw, kw),
                            col_elem.get(kw, None)
                        )
                    if col_elem.get("optional", False):
                        yield (
                            join_kws(prefix, key, col_kw, "enabled"),
                            True
                        )
                    if col_elem.get("colorscale", False):
                        # Use a flag to indicate whether the values
                        # in the column are categorical
                        yield (
                            join_kws(prefix, key, col_kw, "is_categorical"),
                            False
                        )
                        # Use a flag for the color scale to use
                        yield (
                            join_kws(prefix, key, col_kw, "scale"),
                            "Viridis"
                        )
                # If the user is allowed to filter the data
                if elem.get("query", True):
                    for attr in [
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

    def input_selectbox_kwargs(
        self,
        kw,
        options: list,
        names: Optional[list] = None,
        copy_to=None
    ):
        """
        Populate the selectbox element with default kwargs.
        If names is provided, transform the selected option
        into the value in the same index position.
        """

        assert len(options) > 0, f"Empty options list for {kw}"

        if names is not None:
            assert len(options) == len(names), f"Length mismatch for {kw}"
            # Make a dict mapping each name to the values
            names = dict(zip(names, options))

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
            kwargs=dict(
                copy_to=copy_to,
                names=names
            )
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

    def input_value_change(self, kw, names=None, copy_to=None):
        # Get the value provided by the user
        value = st.session_state[self.param_key(kw)]

        # If a name mapping was provided
        if names is not None:

            # The provided value must be a key in the mapping
            assert value in names, f"Invalid value for {kw} ({value})"

            # Update the value to the mapped value
            value = names[value]

        # If this is different than the params
        if value != self.params[kw]:

            # Update the view in the mdata object
            self.update_view_param(kw, value)

            # If a copy_to kwarg was provided
            if copy_to is not None:
                # Modify that value as well
                self.update_view_param(copy_to, value)

    def get_data(self, container: DeltaGenerator):

        # By default, params are complete until proven otherwise
        self.params_complete = True

        if self.params_editable:
            container.write("##### Inputs")

        # Parse the form schema of the object
        self.render_form(container, self.schema)

    def render_form(
        self,
        container: DeltaGenerator,
        schema: Dict[str, dict],
        prefix: str = None
    ):

        # Iterate over the form defined for this view
        for key, elem in schema.items():

            prefix_key = join_kws(prefix, key)

            if elem["type"] == "dataframe":

                self.render_dataframe(prefix, key, elem, container)

            elif elem["type"] == "object":

                if "label" in elem and self.params_editable:
                    container.write(f"**{elem.get('label', prefix_key)}**")

                self.render_form(
                    container,
                    elem["properties"],
                    prefix_key
                )

            elif elem["type"] == "string":

                if self.params_editable:
                    if elem.get("enum") is not None:
                        self.params[prefix_key] = container.selectbox(
                            (
                                key
                                if elem.get("label") is None
                                else elem["label"]
                            ),
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

                if not isinstance(self.params[prefix_key], str):
                    self.params_complete = False

            elif elem["type"] == "float":

                if self.params_editable:
                    self.params[prefix_key] = container.number_input(
                        key if elem.get("label") is None else elem["label"],
                        help=elem.get("help"),
                        **self.input_value_kwargs(prefix_key),
                        **{
                            kw: float(elem[kw])
                            for kw in ["min_value", "max_value"]
                            if elem.get(kw) is not None
                        }
                    )
                if not isinstance(self.params[prefix_key], (float, int)):
                    self.params_complete = False

            elif elem["type"] == "integer":

                # Make sure that the value is an integer
                val = st.session_state.get(
                    self.param_key(prefix_key),
                    elem.get("default")
                )
                if val is None:
                    val = 0

                if not isinstance(val, int):
                    # Update the value
                    self.update_view_param(prefix_key, int(val))

                if self.params_editable:
                    self.params[prefix_key] = container.number_input(
                        key if elem.get("label") is None else elem["label"],
                        help=elem.get("help"),
                        step=1,
                        **self.input_value_kwargs(prefix_key),
                        **{
                            kw: int(elem[kw])
                            for kw in ["min_value", "max_value"]
                            if elem.get(kw) is not None
                        }
                    )

            elif elem["type"] == "boolean":

                if self.params_editable:
                    self.params[prefix_key] = container.checkbox(
                        key if elem.get("label") is None else elem["label"],
                        help=elem.get("help"),
                        **self.input_value_kwargs(prefix_key)
                    )

            elif elem["type"] == "supporting_figure":

                if self.params_editable:
                    # Get the list of all supporting figures
                    all_figures = app.get_supp_figs()

                    if len(all_figures) == 0:
                        container.write("No supporting figures available.")
                        self.params_complete = False

                    else:
                        # Let the user select a supporting figure
                        self.params[prefix_key] = container.selectbox(
                            (
                                key if elem.get("label") is None
                                else elem["label"]
                            ),
                            all_figures,
                            help=elem.get("help"),
                            **self.input_selectbox_kwargs(
                                prefix_key,
                                all_figures
                            )
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

        if not app.has_mdata():
            container.write("No MuData object available.")
            self.params_complete = False
            return

        # Get the orientation (axis) of the data
        axis_kw = join_kws(prefix, key, "axis")

        # Let the user select the orientation to use
        if self.params_editable:
            container.selectbox(
                "Select orientation",
                ["Observations", "Variables"],
                **self.input_selectbox_kwargs(
                    axis_kw,
                    [0, 1],
                    ["Observations", "Variables"]
                )
            )
        axis = self.params[axis_kw]
        assert axis in [0, 1], f"Invalid axis: {axis}"

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
                    axis,
                    container
                )

            df = self.build_dataframe(
                join_kws(prefix, key),
                elem["columns"],
                axis
            )

        # If 'columns' was not specified
        else:

            tables_kw = join_kws(prefix, key, "tables")

            all_tables = app.tree_tables(self.params[axis_kw])

            # If any invalid tables were selected
            if any([
                table not in all_tables
                for table in self.params.get(tables_kw, [])
            ]):
                # Remove any invalid tables
                self.update_view_param(
                    tables_kw,
                    [
                        table
                        for table in self.params.get(tables_kw, [])
                        if table in all_tables
                    ]
                )

            # Let the user select one or more tables
            if self.params_editable:

                container.multiselect(
                    "Select table(s)",
                    all_tables,
                    **self.input_multiselect_kwargs(
                        tables_kw,
                        all_tables
                    )
                )

            # Make a DataFrame with the selected table(s)
            selected_tables: List[str] = self.params.get(tables_kw, [])
            if len(selected_tables) == 0:
                container.write("No tables selected.")
                self.params_complete = False
                return

            df = app.join_dataframe_tables(selected_tables, axis)

        # Let the user optionally filter rows
        if elem.get("query", True):
            filtered_obs = self.render_query(prefix, key, axis, container)
            if filtered_obs is not None:
                df = df.loc[
                    list(set(filtered_obs) & set(df.index))
                ]
                if self.params_editable:
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

            # Let the user select the columns to use
            if self.params_editable:
                # Check to see if the user wants to select all of the columns
                if container.checkbox(
                    "Use all columns",
                    len(selected_columns) == 0
                ):
                    self.update_view_param(selected_columns_kw, [])

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
            if len(selected_columns) > 0:
                df = df[selected_columns].dropna()

        if self.params_editable and df.shape[0] == 0:
            container.write("No data available.")
            return

        self.params[join_kws(prefix, key, "dataframe")] = df
        if self.params_editable:
            container.write(
                "Selected {:,} rows and {:,} columns.".format(*df.shape)
            )

    def render_dataframe_column(
        self,
        prefix: Union[None, str],
        key: str,
        col_kw: str,
        col_elem: dict,
        axis: int,
        container: DeltaGenerator
    ):

        # Set the default values

        # Get the list of tables available for this orientation
        all_tables = app.tree_tables(axis)

        # Table selection
        table_kw = join_kws(prefix, key, col_kw, "table")
        if self.params.get(table_kw) not in all_tables:
            self.update_view_param(
                table_kw,
                app.tree_tables(axis)[0]
            )

        # Get all of the columns for the selected table
        all_cnames = app.list_cnames(
            self.params[table_kw],
            axis=axis
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
        if self.params.get(cname_kw) not in all_cnames:
            cname_val = all_cnames[0]
            self.update_view_param(cname_kw, cname_val)
            self.update_view_param(label_kw, cname_val)

        # If the views are editable
        if self.params_editable:

            # Print the column name
            container.write(f"**{col_elem.get('label', col_kw)}**")

            # Optional column selection
            enabled_kw = join_kws(prefix, key, col_kw, "enabled")
            if col_elem.get("optional", False):
                container.checkbox(
                    "Enabled",
                    **self.input_value_kwargs(enabled_kw)
                )

            # If the column is enabled
            if self.params.get(enabled_kw, True):

                # Make three columns for:
                #   table, column name, and label
                cols = container.columns([1, 1, 1])

                # Select the table of interest
                cols[0].selectbox(
                    "Table",
                    all_tables,
                    **self.input_selectbox_kwargs(table_kw, all_tables)
                )

                # Select the column name
                cols[1].selectbox(
                    "Column",
                    all_cnames,
                    **self.input_selectbox_kwargs(
                        cname_kw,
                        all_cnames,
                        copy_to=label_kw
                    )
                )

                # Input the column label
                cols[2].text_input(
                    "Label",
                    **self.input_value_kwargs(label_kw)
                )

                # If the column is a color column
                if col_elem.get("colorscale", False):

                    # Let the user select whether the colors are categorical
                    is_categorical_kw = join_kws(
                        prefix,
                        key,
                        col_kw,
                        "is_categorical"
                    )

                    # Select a pallete for the colors
                    colorscale_kw = join_kws(prefix, key, col_kw, "scale")

                    container.checkbox(
                        "Categorical Values",
                        help="Color scales may either be categorical or continuous.", # noqa
                        **self.input_value_kwargs(is_categorical_kw)
                    )

                    if self.params[is_categorical_kw]:
                        colors_qualitative = [
                            cname for cname in dir(px.colors.qualitative)
                            if (
                                not cname.startswith("_")
                                and cname != "swatches"
                            )
                        ]

                        container.selectbox(
                            "Select color scale",
                            colors_qualitative,
                            **self.input_selectbox_kwargs(
                                colorscale_kw,
                                colors_qualitative
                            )
                        )

                    else:

                        container.selectbox(
                            "Select color scale",
                            px.colors.named_colorscales(),
                            **self.input_selectbox_kwargs(
                                colorscale_kw,
                                px.colors.named_colorscales()
                            )
                        )

    def build_dataframe(self, key: str, columns: dict, axis: int):
        # Get the information for each column
        return pd.DataFrame({
            col_kw: app.get_dataframe_column(
                axis=axis,
                **{
                    kw: self.param(key, col_kw, kw)
                    for kw in ["table", "cname"]
                }
            )
            for col_kw, col_elem in columns.items()
            if (
                col_elem.get("optional", False) is False
                or
                self.param(key, col_kw, "enabled")
            )
        }).dropna()

    def render_query(
        self,
        prefix: str,
        key: str,
        axis: int,
        container: DeltaGenerator
    ):

        # Set the default values

        # Get the list of tables available for all modalities
        all_tables = app.tree_tables(axis)

        # Table selection
        table_kw = join_kws(prefix, key, "query", "table")
        if self.params.get(table_kw) not in all_tables:
            self.update_view_param(
                table_kw,
                all_tables[0]
            )

        # Get the list of possible columns
        all_cnames = app.list_cnames(
            self.params[table_kw],
            axis=axis
        )

        # Column name selection
        cname_kw = join_kws(prefix, key, "query", "cname")
        if self.params.get(cname_kw) not in all_cnames:

            self.update_view_param(
                cname_kw,
                all_cnames[0]
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
        if self.params_editable:

            # Print the column name
            container.write("#### Filter samples")

            # Make three columns for the table, and column name
            cols = container.columns([1, 1, 1, 1])

            # Select the table of interest
            cols[0].selectbox(
                "Table",
                all_tables,
                **self.input_selectbox_kwargs(table_kw, all_tables)
            )

            # Get the list of possible columns
            all_cnames = app.list_cnames(
                self.params[table_kw],
                axis=axis
            )

            # Select the column name
            cols[1].selectbox(
                "Column",
                all_cnames,
                **self.input_selectbox_kwargs(
                    cname_kw,
                    all_cnames
                )
            )

            # Input the boolean operator
            cols[2].selectbox(
                "Operator",
                [">=", "<=", "==", "!=", ">", "<"],
                **self.input_selectbox_kwargs(
                    expr_kw,
                    [
                        ">=",
                        "<=",
                        "==",
                        "!=",
                        ">",
                        "<"
                    ]
                )
            )

            # Input the boolean value
            cols[3].text_input(
                "Value",
                help="Enter a value to filter samples on.",
                **self.input_value_kwargs(value_kw)
            )

        # Get the values for the query
        query = {
            "table": self.param(key, "query", "table"),
            "cname": self.param(key, "query", "cname"),
            "expr": self.param(key, "query", "expr"),
            "value": self.param(key, "query", "value")
        }

        # If no value is provided
        if query['value'] is None or len(query['value']) == 0:
            if self.params_editable:
                container.write("Provide a value to filter samples.")
            return

        if query["table"] == "Observation Metadata":
            query["table"] = "obs"
            query["modality"] = None

        else:
            query["modality"], query["table"] = query["table"].split(".")

        # Get the table
        table = app.get_dataframe_table(
            query["modality"],
            query["table"],
            axis
        )
        if table is None:
            if self.params_editable:
                container.write("No data available for filtering.")
            return

        msg = f"Column {query['cname']} not found in table."
        assert query['cname'] in table.columns, msg

        # Apply the filter
        try:
            table = table.query("{cname} {expr} {value}".format(**query))
        except Exception as e:
            container.write("Error while filtering")
            container.exception(e)
            return

        # If no values are returned
        if table.shape[0] == 0:
            container.write("No samples match the filter criteria.")
            return

        return table.index
