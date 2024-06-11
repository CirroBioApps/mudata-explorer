from mudata_explorer import app
from mudata_explorer.helpers.join_kws import join_kws
from mudata_explorer.base import all_transforms, get_transform
from typing import Dict, List, Optional
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
    # Optional help text
    help_text: Optional[str] = False

    def param(self, *kws, default=None):
        return self.params.get(join_kws(*kws), default)

    def get_schema_defaults(self, schema: dict, prefix=None):
        """Yield the default values for each item in the schema."""

        for elem_key, elem in schema.items():
            assert isinstance(elem, dict), f"Expected dict, got {type(elem)}"
            assert "type" in elem, f"Missing 'type' key in schema: {elem}"

            # Join the prefix and key
            key = join_kws(prefix, elem_key)

            if elem["type"] == "object":
                yield from self.get_schema_defaults(elem["properties"], key)

            elif elem["type"] in [
                "string",
                "float",
                "boolean",
                "integer",
                "supporting_figure"
            ]:
                yield (
                    key,
                    elem.get("default", None)
                )
                if elem.get("optional", False):
                    yield (
                        join_kws(key, "enabled"),
                        True
                    )

            elif elem["type"] == "dataframe":
                # Each dataframe may be oriented to the obs or var
                yield (
                    join_kws(key, "axis"),
                    elem.get("axis", 0)
                )
                # If column selection is not enabled,
                # the user will select >= 1 tables
                if len(elem.get("columns", {})) == 0:
                    yield (
                        join_kws(key, "tables"),
                        elem.get("tables", [])
                    )

                # Add any columns specified in the schema
                for col_kw, col_elem in elem.get("columns", {}).items():
                    for kw in ["table", "cname", "label"]:
                        yield (
                            join_kws(key, col_kw, kw),
                            col_elem.get(kw, None)
                        )
                    if col_elem.get("optional", False):
                        yield (
                            join_kws(key, col_kw, "enabled"),
                            True
                        )
                    if col_elem.get("colorscale", False):
                        # Use a flag to indicate whether the values
                        # in the column are categorical
                        yield (
                            join_kws(key, col_kw, "is_categorical"),
                            col_elem.get("is_categorical", False)
                        )
                        # Use a flag for the color scale to use
                        yield (
                            join_kws(key, col_kw, "scale"),
                            "D3" if col_elem.get("is_categorical", False) else "Viridis"
                        )

                # Filtering of the rows and columns
                for axis_kw in ["rows_query", "cols_query"]:
                    for attr in [
                        "type",   # 'value' or 'index'
                        "table",  # Name of the table for value filtering
                        "cname",  # Name of the column for value filtering
                        "expr",   # Type of comparison ('in', 'not in', '==', '!=', '>', '<', '>=', '<=') # noqa
                        "value"   # Either the value or the specific indices
                    ]:
                        yield (join_kws(key, axis_kw, "query", attr), "")

                # Transforming the values
                yield (
                    join_kws(key, "transforms"),
                    elem.get("transforms", [])
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
        if self.ix == -1:
            mdata.uns["mudata-explorer-process"]["params"][kw] = value
        else:
            mdata.uns["mudata-explorer-views"][self.ix]["params"][kw] = value

        # Save the MuData object
        app.set_mdata(mdata)

        # Also update the params object
        self.params[kw] = value

    def delete_view_param(self, kw):
        # Get the MuData object
        mdata = app.get_mdata()

        # Delete the value of this param for this view
        if self.ix == -1:
            params = mdata.uns["mudata-explorer-process"]["params"]
        else:
            params = mdata.uns["mudata-explorer-views"][self.ix]["params"]
        if kw in params:
            del params[kw]

        # Save the MuData object
        app.set_mdata(mdata)

        # Also update the params object
        if kw in self.params:
            del self.params[kw]

    def input_selectbox_kwargs(
        self,
        kw,
        options: list,
        names: Optional[list] = None,
        copy_to=None,
        invalidate=[]
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
                names=names,
                invalidate=invalidate
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

    def input_value_change(
        self,
        kw,
        names=None,
        copy_to=None,
        invalidate=[]
    ):
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

            for kw in invalidate:
                self.delete_view_param(kw)

    def get_data(self, container: Optional[DeltaGenerator] = None):

        if container is None:
            container = st.container()

        # By default, params are complete until proven otherwise
        self.params_complete = True

        if self.params_editable:
            # Show the name of the view
            container.write(
                f"#### {self.ix + 1}. {self.name}"
                if self.ix != -1
                else f"#### {self.name}"
            )

        # Parse the form schema of the object
        self.render_form(container, self.schema)

        if self.params_editable and self.ix >= 0:
            container.button(
                "Save Changes",
                key=f"save-changes-{self.ix}",
                on_click=self.save_changes
            )

    def save_changes(self):
        del st.query_params["edit-view"]

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

                self.render_dataframe(
                    prefix_key,
                    elem,
                    container.container(border=True)
                )

            elif elem["type"] == "object":

                if "label" in elem and self.params_editable:
                    container.write(f"**{elem.get('label', prefix_key)}**")

                self.render_form(
                    container,
                    elem["properties"],
                    prefix_key
                )
                if self.params_editable:
                    container.write("---")

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
        key: str,
        elem: dict,
        container: DeltaGenerator
    ):

        if "label" in elem and self.params_editable:
            container.write(f"**{elem.get('label')}**")

        if not app.has_mdata():
            container.write("No MuData object available.")
            self.params_complete = False
            return

        # Get the orientation (axis) of the data
        axis_kw = join_kws(key, "axis")

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

        # Tables are defined either by selecting specific columns,
        # or by selecting one or more tables

        # If 'columns' were specified
        if len(elem.get("columns", {})) > 0:

            # Build the DataFrame from the specified columns
            df = self.render_dataframe_columns(
                key,
                elem["columns"],
                axis,
                container
            )

        # If 'columns' was not specified
        else:

            # Build the DataFrame from the selected tables
            df = self.render_dataframe_tables(
                key,
                axis,
                container
            )

            if df is not None:

                # The user can filter the data along the columns
                df = self.filter_dataframe_cols(key, axis, df, container)

        if df is not None:

            # The user can filter the data along the rows
            df = self.filter_dataframe_rows(key, axis, df, container)

            # The user can transform the values in the DataFrame
            df = self.transform_dataframe(key, df, container)

            # Drop any columns which are entirely missing
            dropped_cols = df.shape[1] - df.dropna(axis=1, how="all").shape[1]
            df = df.dropna(axis=1, how="all")
            if dropped_cols:
                if self.params_editable:
                    container.write(
                        f"Removed {dropped_cols:,} columns with entirely missing values." # noqa
                    )

            # Drop any null values
            dropped_rows = df.shape[0] - df.dropna().shape[0]
            df = df.dropna()
            if dropped_rows:
                if self.params_editable:
                    container.write(
                        f"Removed {dropped_rows:,} rows with missing values."
                    )

        if self.params_editable and (df is None or df.shape[0] == 0):
            container.write("No data available.")
            return

        self.params[join_kws(key, "dataframe")] = df
        if self.params_editable and df is not None:
            container.write(
                "Selected {:,} rows and {:,} columns.".format(*df.shape)
            )

    def render_dataframe_columns(
        self,
        key: str,
        columns: dict,
        axis: int,
        container: DeltaGenerator
    ):
        # Iterate over the columns in the schema and prompt the user
        # for which data to supply for each
        for col_kw, col_elem in columns.items():

            self.render_dataframe_column(
                key,
                col_kw,
                col_elem,
                axis,
                container
            )

        return self.build_dataframe(
            key,
            columns,
            axis
        )

    def render_dataframe_column(
        self,
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
        table_kw = join_kws(key, col_kw, "table")
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
        label_kw = join_kws(key, col_kw, "label")
        if self.params.get(label_kw) is None:
            self.update_view_param(
                label_kw,
                "Label"
            )

        # Column name selection
        cname_kw = join_kws(key, col_kw, "cname")
        if self.params.get(cname_kw) not in all_cnames:
            cname_val = all_cnames[0]
            self.update_view_param(cname_kw, cname_val)
            self.update_view_param(label_kw, cname_val)

        # If the views are editable
        if self.params_editable:

            # Print the column name
            container.write(f"**{col_elem.get('label', col_kw)}**")

            # Optional column selection
            enabled_kw = join_kws(key, col_kw, "enabled")
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
                        key,
                        col_kw,
                        "is_categorical"
                    )

                    # Select a pallete for the colors
                    colorscale_kw = join_kws(key, col_kw, "scale")

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

    def render_dataframe_tables(
        self,
        key,
        axis: int,
        container: DeltaGenerator
    ):

        tables_kw = join_kws(key, "tables")

        all_tables = app.tree_tables(axis)

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

        return app.join_dataframe_tables(selected_tables, axis)

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

    def filter_dataframe_cols(
        self,
        key: str,
        axis: int,
        df: pd.DataFrame,
        parent_container: DeltaGenerator
    ) -> pd.DataFrame:
        """
        Let the user select a subset of columns for analysis.
        """
        if self.ix == -1:
            container = parent_container.expander("Filter Columns")
        else:
            container = parent_container.container(border=True)
            container.write("**Filter Columns**")
        df = self.render_query(
            join_kws(key, "cols_query"),
            0 if axis else 1,  # Filter the opposite axis
            1,  # Perform the filtering on the columns
            df,
            container
        )
        if self.params_editable and df is not None:
            container.write(f"Number of filtered columns: {df.shape[1]:,}")
        return df

    def filter_dataframe_rows(
        self,
        key: str,
        axis: int,
        df: pd.DataFrame,
        parent_container: DeltaGenerator
    ) -> pd.DataFrame:
        """
        Let the user select a subset of rows for analysis.
        """
        if self.ix == -1:
            container = parent_container.expander("Filter Rows")
        elif self.params_editable:
            container = parent_container.container(border=True)
            container.write("**Filter Rows**")
        else:
            container = parent_container

        df = self.render_query(
            join_kws(key, "rows_query"),
            axis,
            0,  # Perform the filtering on the rows
            df,
            container
        )
        if self.params_editable and df is not None:
            container.write(f"Number of filtered rows: {df.shape[0]:,}")
        return df

    def render_query(
        self,
        key: str,
        axis: int,
        filter_axis: int,
        df: pd.DataFrame,
        container: DeltaGenerator
    ) -> pd.DataFrame:

        # Set the default values

        # Type of selection, either by value or by index
        type_kw = join_kws(key, "query", "type")
        if self.params.get(type_kw) is None:
            # By default, select by value
            self.update_view_param(
                type_kw,
                "value"
            )

        # Get the list of tables available for all modalities
        all_tables = app.tree_tables(axis)

        # Table selection
        table_kw = join_kws(key, "query", "table")
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
        cname_kw = join_kws(key, "query", "cname")
        if self.params.get(cname_kw) not in all_cnames:

            self.update_view_param(
                cname_kw,
                all_cnames[0]
            )

        # Boolean operator selection
        expr_kw = join_kws(key, "query", "expr")
        if self.params.get(expr_kw) is None:

            self.update_view_param(
                expr_kw,
                ">="
            )

        # Comparison value selection
        # There are a few cases where we will set a new default value
        value_kw = join_kws(key, "query", "value")
        if (
            # if no value is set
            self.params.get(value_kw) is None
            or (
                # If the value is a list but the
                # expression is not 'in' or 'not in'
                isinstance(self.params.get(value_kw), list)
                and self.params.get(expr_kw) not in ["in", "not in"]
            )
            or (
                # If the value is not a list but the
                # expression is 'in' or 'not in'
                not isinstance(self.params.get(value_kw), list)
                and self.params.get(expr_kw) in ["in", "not in"]
            )
        ):

            self.update_view_param(
                value_kw,
                (
                    []
                    if self.params.get(expr_kw) in ["in", "not in"]
                    else ""
                )
            )

        # If the views are editable
        if self.params_editable:

            # First figure out whether we're selecting by value
            # or by selecting specific indices
            names = [
                "Filtering by Value",
                "Selecting Specific {}".format(
                    "Rows" if filter_axis == 0 else "Columns"
                )
            ]
            values = ["value", "index"]
            container.selectbox(
                "Select by:",
                names,
                **self.input_selectbox_kwargs(
                    type_kw,
                    values,
                    names,
                    invalidate=[value_kw]
                )
            )

            # If we're selecting by value
            if self.params[type_kw] == "value":

                # Select the table of interest
                container.selectbox(
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
                container.selectbox(
                    "Column",
                    all_cnames,
                    **self.input_selectbox_kwargs(
                        cname_kw,
                        all_cnames
                    )
                )
                operators = [">=", "<=", "==", "!=", ">", "<", "in", "not in"]

            else:
                operators = ["in", "not in"]

            # Input the comparison operator
            container.selectbox(
                "Keep items which are:",
                operators,
                **self.input_selectbox_kwargs(expr_kw, operators)
            )

            if self.params[expr_kw] in ["in", "not in"]:

                # If we're selecting by value
                if self.params[type_kw] == "value":

                    # Let the user select specific values, selecting
                    # from the values in the table
                    value_options = app.get_dataframe_column(
                        axis,
                        self.params[table_kw],
                        self.params[cname_kw]
                    ).unique()

                # Otherwise, we're selecting by index
                else:

                    # Get the index values along the axis
                    # which will be filtered
                    value_options = (
                        df.index.values if filter_axis == 0
                        else df.columns.values
                    )

                # Input the comparison value directly
                container.multiselect(
                    "Value",
                    value_options,
                    help="Select the values for filtering",
                    **self.input_multiselect_kwargs(
                        value_kw,
                        value_options
                    )
                )

            else:

                # Input the comparison value directly
                container.text_input(
                    "Value",
                    help="Enter a value to filter samples on.",
                    **self.input_value_kwargs(value_kw)
                )

        # Get the values for the query
        query = {
            kw: self.param(key, "query", kw)
            for kw in [
                "type",
                "table",
                "cname",
                "expr",
                "value"
            ]
        }

        # If the user is selecting certain indices
        if query["type"] == "index":

            # Only keep the values which are in that index
            query["value"] = [
                val for val in query["value"]
                if val in (
                    df.index if filter_axis == 0
                    else df.columns
                )
            ]

        # If no value is provided
        if query['value'] is None or len(query['value']) == 0:
            if self.params_editable:
                container.write("Provide a value to filter samples.")
            return df

        # If the user is selecting certain indices
        if query["type"] == "index":

            # If we are dropping rows/cols
            assert query["expr"] in ["in", "not in"]
            if query["expr"] == "not in":
                return df.drop(
                    query["value"],
                    axis=filter_axis
                )
            else:
                return df.reindex(
                    query["value"],
                    axis=filter_axis
                )

        # If we are filtering by value

        if query["table"] == "Observation Metadata":
            query["table"] = "metadata"
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
            return df

        msg = f"Column {query['cname']} not found in table."
        assert query['cname'] in table.columns, msg

        # If we are selecting certain values
        if query["expr"] in ["in", "not in"]:
            query["value"] = [
                val for val in query["value"]
                if val in table[query["cname"]].values
            ]
            if len(query["value"]) == 0:
                if self.params_editable:
                    container.write("No values match the filter criteria.")
                return df

            if query["expr"] == "in":
                return df.loc[
                    df[query["cname"]].isin(query["value"])
                ]
            else:
                return df.loc[
                    ~df[query["cname"]].isin(query["value"])
                ]

        # Apply the filter, trying both string and numeric values
        filtered_table = self._filter_table(table, query)

        # If no values are returned
        if filtered_table is None:
            container.write("No samples match the filter criteria.")
            return df

        # Subset the larger table by the filtered table
        if filter_axis == 0:
            return df.loc[list(set(filtered_table.index) & set(df.index))]
        else:
            return df[list(set(filtered_table.index) & set(df.columns))]

    @staticmethod
    def _filter_table(table, query):

        filtered_table = None

        try:
            filtered_table = table.query(
                "{cname} {expr} {value}".format(**query)
            )
        except: # noqa
            pass

        if isinstance(filtered_table, pd.DataFrame):
            if filtered_table.shape[0] > 0:
                return filtered_table

        try:
            filtered_table = table.query(
                "{cname} {expr} '{value}'".format(**query)
            )
        except: # noqa
            pass

        if isinstance(filtered_table, pd.DataFrame):
            if filtered_table.shape[0] > 0:
                return filtered_table

    def _get_values_in_column(
        self,
        axis: int,
        table: str,
        cname: str
    ) -> List[str]:
        """
        Get the unique values in the specified column.
        """
        # Get the unique values
        return app.get_dataframe_column(axis, table, cname).unique()

    def transform_dataframe(
        self,
        key: str,
        df: pd.DataFrame,
        parent_container: DeltaGenerator
    ) -> pd.DataFrame:

        if self.ix == -1:
            container = parent_container.expander("Transform Values")
        elif self.params_editable:
            container = parent_container.container(border=True)
            container.write("**Transform Values**")
        else:
            container = parent_container

        # Get the list of transformations
        transforms = self.param(key, "transforms")

        if self.params_editable:

            # Show each of the transformations, and also allow
            # the user to remove those transformations
            for ix, transform in enumerate(transforms):
                cols = container.columns([4, 1])
                cols[0].button(
                    f"{ix + 1}: {get_transform(transform).name}",
                    use_container_width=True
                )
                cols[1].button(
                    "Remove",
                    key=f"remove_transform_{ix}",
                    on_click=self._remove_transform,
                    args=(key, ix),
                    use_container_width=True
                )

            # Let the user add a transformation
            if container.button(
                "Add Transformation",
                key=f"add-transformation-{key}-{self.ix}"
            ):
                self._add_transform_popup(key)

        # Run all of the transformations
        for transform in transforms:
            try:
                df = get_transform(transform).run(df)
            except Exception as e:
                container.exception(e)

        return df

    def _remove_transform(self, key, ix):
        """Remove a particular transform from the list."""

        # Get the list of transformations
        transforms: list = self.param(key, "transforms")

        # Remove the element
        transforms.pop(ix)

        # Update the list of transformations
        self.update_view_param(
            join_kws(key, "transforms"),
            transforms
        )

    @st.experimental_dialog("Add Transformation")
    def _add_transform_popup(self, key):
        st.selectbox(
            "Select Transformation",
            [""] + [
                transform.name
                for transform in all_transforms().values()
            ],
            index=0,
            key=f"add_transform_{key}",
        )
        if st.session_state[f"add_transform_{key}"] != "":
            if st.button("Add"):
                self._add_transform(key)
                st.rerun()

    def _add_transform(self, key):
        """Add a transform to the list."""

        # Get the selected transform (by name)
        selected_name = st.session_state[f"add_transform_{key}"]

        # Stop if no name was selected
        if len(selected_name) == 0:
            return

        # Get the id of that transform
        selected_id = None
        for transform_id, transform in all_transforms().items():
            if transform.name == selected_name:
                selected_id = transform_id
                break

        msg = f"Couldn't find transform: {selected_name}"
        assert selected_id is not None, msg

        # Get the list of transformations
        transforms: list = self.param(key, "transforms")

        # Add the new transform
        transforms.append(selected_id)

        # Update the list of transformations
        self.update_view_param(
            join_kws(key, "transforms"),
            transforms
        )
