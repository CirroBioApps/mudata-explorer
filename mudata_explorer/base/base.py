import json
from mudata_explorer.helpers.join_kws import join_kws
from mudata_explorer.base import all_transforms, get_transform
from mudata_explorer.app.mdata import get_mdata, set_mdata, has_mdata
from mudata_explorer.app.mdata import get_supp_figs
from mudata_explorer.app.mdata import tree_tables, list_cnames, join_dataframe_tables, get_dataframe_column
from mudata_explorer.helpers.views import get_views
from mudata_explorer.app.process import nest_params
from mudata_explorer.app.query_params import get_edit_views_flag
from typing import Dict, List, Optional
import pandas as pd
import plotly.express as px
import streamlit as st
from streamlit.delta_generator import DeltaGenerator


class MuDataAppHelpers:

    ix: int
    params: dict = {}
    uns: dict = {}
    schema: dict
    type: str
    name: str
    category: str
    # Flag used to indicate whether all parameters have been provided
    params_complete: bool
    # Flag used to indicate whether the parameters are editable
    params_editable: bool
    # Optional help text
    help_text: Optional[str] = False
    # Optionally attach a MuData object to the object
    mdata = None

    def param(self, *kws, default=None):
        return self.params.get(join_kws(*kws), default)

    def get_schema_defaults(self, schema: dict, prefix=None):
        """Yield the default values for each item in the schema."""

        for elem_key, elem in schema.items():
            assert isinstance(elem, dict), f"Expected dict, got {type(elem)}"
            assert "type" in elem, f"Missing 'type' key in schema: {elem}"

            # Join the prefix and key
            key = join_kws(prefix, elem_key)

            # By default, nothing is in the sidebar
            yield (
                join_kws(key, "sidebar"),
                False
            )

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
                    yield (
                        join_kws(key, "enabled", "sidebar"),
                        False
                    )

            elif elem["type"] == "dataframe":

                if elem.get("optional", False):
                    yield (
                        join_kws(key, "enabled"),
                        True
                    )
                    yield (
                        join_kws(key, "enabled", "sidebar"),
                        False
                    )

                # Each dataframe may be oriented to the obs or var
                yield (
                    join_kws(key, "axis"),
                    elem.get("axis", 0)
                )
                yield (
                    join_kws(key, "axis", "sidebar"),
                    False
                )

                # If column selection is not enabled,
                # the user will select >= 1 tables
                if len(elem.get("columns", {})) == 0:
                    yield (
                        join_kws(key, "tables"),
                        elem.get("tables", [])
                    )
                    yield (
                        join_kws(key, "tables", "sidebar"),
                        False
                    )

                # Add any columns specified in the schema
                for col_kw, col_elem in elem.get("columns", {}).items():
                    for kw in ["table", "cname", "label"]:
                        yield (
                            join_kws(key, col_kw, kw),
                            col_elem.get(kw, None)
                        )
                        yield (
                            join_kws(key, col_kw, kw, "sidebar"),
                            False
                        )

                    if col_elem.get("optional", False):
                        yield (
                            join_kws(key, col_kw, "enabled"),
                            True
                        )
                        yield (
                            join_kws(key, col_kw, "enabled", "sidebar"),
                            False
                        )

                    if col_elem.get("colorscale", False):
                        # Use a flag to indicate whether the values
                        # in the column are categorical
                        yield (
                            join_kws(key, col_kw, "is_categorical"),
                            col_elem.get("is_categorical", False)
                        )
                        yield (
                            join_kws(key, col_kw, "is_categorical", "sidebar"),
                            False
                        )

                        # Use a flag for the color scale to use
                        yield (
                            join_kws(key, col_kw, "scale"),
                            (
                                "D3"
                                if col_elem.get("is_categorical", False)
                                else "Viridis"
                            )
                        )
                        yield (
                            join_kws(key, col_kw, "scale", "sidebar"),
                            False
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
                        yield (
                            join_kws(key, axis_kw, "query", "sidebar"),
                            False
                        )

                # Transforming the values
                yield (
                    join_kws(key, "transforms"),
                    elem.get("transforms", [])
                )
                yield (
                    join_kws(key, "transforms", "sidebar"),
                    False
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
        # Get the element in the global mdata object for this view
        mdata, elem = self._mdata_elem()

        # Modify the value of this param for this view
        elem["params"][kw] = value

        # Save the MuData object
        set_mdata(mdata)

        # Also update the params object
        self.params[kw] = value

    def _mdata_elem(self):
        """Return the element in the MuData object corresponding to this view."""
        # Get the MuData object
        mdata = get_mdata()

        # Point to the element to modify
        if self.ix == -1:
            elem = mdata.uns["mudata-explorer-process"]
        else:
            elem = mdata.uns["mudata-explorer-views"][self.ix]
        
        return mdata, elem
    
    def delete_view_param(self, kw):
        # Get the element in the global mdata object for this view
        mdata, elem = self._mdata_elem()

        # Delete the value of this param for this view
        if kw in elem["params"]:
            del elem["params"][kw]

            # Save the MuData object
            set_mdata(mdata)

        # Also update the params object
        if kw in self.params:
            del self.params[kw]

    def update_view_uns(self, kw, value):

        # Get the element in the global mdata object for this view
        mdata, elem = self._mdata_elem()

        # Make sure that the uns key exists
        if "uns" not in elem:
            elem["uns"] = {}

        # Modify the value of this uns for this view
        elem["uns"][kw] = value

        # Save the MuData object
        set_mdata(mdata)

        # Also update the params object
        self.uns[kw] = value

    def delete_view_uns(self, kw):
        # Get the element in the global mdata object for this view
        mdata, elem = self._mdata_elem()

        # Make sure that the uns key exists
        if "uns" not in elem:
            elem["uns"] = {}

        if kw in elem["uns"]:
            del elem["uns"][kw]

            # Save the MuData object
            set_mdata(mdata)

        # Also update the uns object
        if kw in self.uns:
            del self.uns[kw]

    def _input_container(self, kw: str) -> DeltaGenerator:
        """
        Return a container which can be used for an input element.
        If we're in the menu where all params are editable, then
        a checkbox will be displayed to the side which will toggle
        whether the param will be displayed in the sidebar.
        """
        # If we're in the menu where all params are editable
        if self.params_editable:
            # Set up two columns, one for the input and one to toggle
            # whether this param is viewable in the sidebar
            input, toggle = st.columns([2, 1])
        else:
            input = st.container()

        if self.params_editable:
            with toggle:
                self.toggle_sidebar(kw)

        return input

    def selectbox(
        self,
        label: str,
        options: list,
        kw: str,
        display_options: Optional[list]=None,
        copy_to=None,
        invalidate=[],
        help=None
    ):
        display_options = display_options if display_options else options

        return self._input_container(kw).selectbox(
            label,
            display_options,
            help=help,
            **self.input_selectbox_kwargs(
                kw,
                options,
                display_options,
                copy_to=copy_to,
                invalidate=invalidate
            )
        )

    def multiselect(
        self,
        label: str,
        options: list,
        kw: str,
        copy_to=None,
        help=None
    ):

        return self._input_container(kw).multiselect(
            label,
            options=options,
            help=help,
            **self.input_multiselect_kwargs(
                kw,
                options,
                copy_to=copy_to
            )
        )

    def checkbox(self, label, kw, **kwargs):

        return self._input_container(kw).checkbox(label, **kwargs)

    def number_input(
        self,
        label: str,
        kw: str,
        help=None,
        copy_to=None,
        **kwargs
    ):

        return self._input_container(kw).number_input(
            label,
            help=help,
            **self.input_value_kwargs(kw, copy_to=copy_to),
            **kwargs
        )

    def text_input(self, label: str, kw: str, help=None, copy_to=None):

        return self._input_container(kw).text_input(
            label,
            help=help,
            **self.input_value_kwargs(kw, copy_to=copy_to)
        )

    def text_area(self, label: str, kw: str, help=None, copy_to=None):

        return self._input_container(kw).text_area(
            label,
            help=help,
            **self.input_value_kwargs(kw, copy_to=copy_to)
        )

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

    def get_data(self):

        # By default, params are complete until proven otherwise
        self.params_complete = True

        if self.params_editable:
            # Show the name of the view
            st.write(
                f"#### {self.ix + 1}. {self.name}"
                if self.ix != -1
                else f"#### {self.name}"
            )
            if st.query_params.get("edit-view"):
                st.write(self.help_text)

        # Parse the form schema of the object
        self.render_form(self.schema)

        if self.params_editable and self.ix >= 0:

            # Now make the display, catching any errors
            try:
                self.display()
            except Exception as e:
                # Log the full traceback of the exception
                st.exception(e)

            st.button(
                ":information_source: SDK Snippet",
                on_click=_show_view_sdk_snippet,
                key=f"show-view-sdk-snippet-{self.ix}",
                args=(self.ix,),
                help="Show an example for configuration via SDK."
            )

            st.button(
                ":page_facing_up: Save Changes",
                key=f"save-changes-{self.ix}",
                on_click=self.save_changes
            )

    def save_changes(self):
        del st.query_params["edit-view"]

    def render_form(
        self,
        schema: Dict[str, dict],
        prefix: str = None
    ):

        # Iterate over the form defined for this view
        for key, elem in schema.items():

            prefix_key = join_kws(prefix, key)

            # The element will be displayed if the global
            # param self.params_editiable is True (indicating that
            # the full set of parameters will be displayed)

            if elem["type"] == "dataframe":

                self.render_dataframe(prefix_key, elem)

            elif elem["type"] == "object":

                if "label" in elem and self.show_param(prefix_key):
                    st.write(f"**{elem.get('label', prefix_key)}**")

                if elem.get("help") is not None and self.show_param(prefix_key):
                    st.write(elem.get('help'))

                self.render_form(
                    elem["properties"],
                    prefix_key
                )
                if self.show_param(prefix_key):
                    st.write("---")

            elif elem["type"] == "string":

                if self.show_param(prefix_key):
                    if elem.get("enum") is not None:
                        self.params[prefix_key] = self.selectbox(
                            label=(
                                key
                                if elem.get("label") is None
                                else elem["label"]
                            ),
                            options=elem["enum"],
                            help=elem.get("help"),
                            kw=prefix_key
                        )

                    elif elem.get("multiline", False):
                        self.params[prefix_key] = self.text_area(
                            label=(
                                key
                                if elem.get("label") is None
                                else elem["label"]
                            ),
                            kw=prefix_key,
                            help=elem.get("help")
                        )

                    else:
                        self.params[prefix_key] = self.text_input(
                            label=(
                                key
                                if elem.get("label") is None
                                else elem["label"]
                            ),
                            kw=prefix_key,
                            help=elem.get("help")
                        )

                if not isinstance(self.params[prefix_key], str):
                    st.write(f"Please provide a valid string ({prefix_key}).")
                    self.params_complete = False

            elif elem["type"] == "float":

                if self.show_param(prefix_key):
                    self.params[prefix_key] = self.number_input(
                        label=key if elem.get("label") is None else elem["label"],
                        help=elem.get("help"),
                        kw=prefix_key,
                        **{
                            kw: float(elem[kw])
                            for kw in ["min_value", "max_value", "step"]
                            if elem.get(kw) is not None
                        }
                    )
                if not isinstance(self.params[prefix_key], (float, int)):
                    self.params_complete = False

            elif elem["type"] == "integer":

                if self.show_param(prefix_key):
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

                    self.params[prefix_key] = self.number_input(
                        label=key if elem.get("label") is None else elem["label"],
                        kw=prefix_key,
                        help=elem.get("help"),
                        step=1,
                        **{
                            kw: int(elem[kw])
                            for kw in ["min_value", "max_value"]
                            if elem.get(kw) is not None
                        }
                    )

            elif elem["type"] == "boolean":

                if self.show_param(prefix_key):
                    self.params[prefix_key] = self.checkbox(
                        key if elem.get("label") is None else elem["label"],
                        kw=prefix_key,
                        help=elem.get("help"),
                        **self.input_value_kwargs(prefix_key)
                    )

            elif elem["type"] == "supporting_figure":

                if self.show_param(prefix_key):
                    # Get the list of all supporting figures
                    all_figures = get_supp_figs()

                    if len(all_figures) == 0:
                        st.write("No supporting figures available.")
                        self.params_complete = False

                    else:
                        # Let the user select a supporting figure
                        self.params[prefix_key] = self.selectbox(
                            label=(
                                key if elem.get("label") is None
                                else elem["label"]
                            ),
                            options=all_figures,
                            help=elem.get("help"),
                            kw=prefix_key
                        )

            else:
                raise Exception(f"Unsupported type: {elem['type']}")

    def render_dataframe(
        self,
        key: str,
        elem: dict
    ):

        if "label" in elem and self.show_param(key):
            st.write(f"**{elem.get('label')}**")

        if elem.get("help") is not None and self.show_param(key):
            st.write(elem.get('help'))

        # Optional column selection
        enabled_kw = join_kws(key, "enabled")
        if elem.get("optional", False) and self.show_param(enabled_kw):
            self.checkbox(
                "Enabled",
                kw=enabled_kw,
                **self.input_value_kwargs(enabled_kw)
            )

        # If the column is disabled
        if not self.params.get(enabled_kw, True):
            return

        if self.show_param(enabled_kw) and has_mdata() is False:
            st.write("No MuData object available.")
            self.params_complete = False
            return

        # Get the orientation (axis) of the data
        axis_kw = join_kws(key, "axis")

        # Let the user select the orientation to use
        if self.show_param(axis_kw):
            self.selectbox(
                label="Select orientation",
                options=[0, 1],
                kw=axis_kw,
                display_options=["Observations", "Variables"]
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
                axis
            )

        # If 'columns' was not specified
        else:

            # Build the DataFrame from the selected tables
            df = self.render_dataframe_tables(
                key,
                axis
            )

            if df is not None:

                # The user can filter the data along the columns
                df = self.filter_dataframe_cols(key, axis, df)

        if df is not None:

            # The user can filter the data along the rows
            df = self.filter_dataframe_rows(key, axis, df)

            # The user can transform the values in the DataFrame
            df = self.transform_dataframe(key, df)

            # Drop any columns which are entirely missing
            dropped_cols = df.shape[1] - df.dropna(axis=1, how="all").shape[1]
            df = df.dropna(axis=1, how="all")
            if dropped_cols:
                if self.show_param(key):
                    st.write(
                        f"Removed {dropped_cols:,} columns with entirely missing values." # noqa
                    )

            if elem.get("dropna", True):

                # Drop any null values
                n = df.shape[0] - df.dropna().shape[0]
                df = df.dropna()
                if n:
                    if self.show_param(key):
                        msg = f"Removed {n:,} rows with missing values."
                        st.write(msg)

        if self.show_param(key) and (df is None or df.shape[0] == 0):
            st.write("No data available.")
            return

        self.params[join_kws(key, "dataframe")] = df
        if self.show_param(key) and df is not None:
            st.write(
                "Selected {:,} rows and {:,} columns.".format(*df.shape)
            )

    def render_dataframe_columns(
        self,
        key: str,
        columns: dict,
        axis: int
    ):
        # Iterate over the columns in the schema and prompt the user
        # for which data to supply for each
        for col_kw, col_elem in columns.items():

            self.render_dataframe_column(
                join_kws(key, col_kw),
                col_elem,
                axis
            )

        return self.build_dataframe(
            key,
            columns,
            axis
        )

    def render_dataframe_column(
        self,
        col_kw: str,
        col_elem: dict,
        axis: int
    ):
        # Set up the kw elements for this column
        table_kw = join_kws(col_kw, "table")
        label_kw = join_kws(col_kw, "label")
        cname_kw = join_kws(col_kw, "cname")
        enabled_kw = join_kws(col_kw, "enabled")
        is_categorical_kw = join_kws(col_kw, "is_categorical")
        colorscale_kw = join_kws(col_kw, "scale")

        # Set the default values
        if self.mdata is None:

            # Get the list of tables available for this orientation
            all_tables = tree_tables(axis)

            # Table selection
            # Only select tables which are valid for this axis
            if self.params.get(table_kw) is None:
                self.update_view_param(table_kw, [])
            if any([
                table not in all_tables
                for table in self.params.get(table_kw)
            ]):
                self.update_view_param(
                    table_kw,
                    [
                        table
                        for table in self.params.get(table_kw)
                        if table in all_tables
                    ]
                )

            # Get all of the columns for the selected table
            all_cnames = list_cnames(
                self.params[table_kw],
                axis=axis
            )

            # Column label selection
            if self.params.get(label_kw) is None:
                self.update_view_param(
                    label_kw,
                    "Label"
                )

            # Column name selection
            if self.params.get(cname_kw) not in all_cnames:
                cname_val = (
                    all_cnames[0]
                    if len(all_cnames) > 0
                    else ""
                )
                self.update_view_param(cname_kw, cname_val)
                self.update_view_param(label_kw, cname_val)

        # If any of the input elements for this selector are viewable
        if self.show_param(col_kw, enabled_kw, table_kw, label_kw, cname_kw):

            # Print the column name
            st.write(f"**{col_elem.get('label', col_kw)}**")

        # Optional column selection
        if col_elem.get("optional", False):
            if self.show_param(enabled_kw):
                self.checkbox(
                    "Enabled",
                    kw=enabled_kw,
                    **self.input_value_kwargs(enabled_kw)
                )

        # If the column is enabled
        if self.params.get(enabled_kw, True):

            # Select the table of interest
            if self.show_param(table_kw):
                self.multiselect(
                    label="Table",
                    options=all_tables,
                    kw=table_kw,
                    help="Select the table(s) containing input data"
                )

                # If no tables have been selected
                if len(self.params.get(table_kw, [])) == 0:
                    st.write("No tables selected.")
                    self.params_complete = False
                    return

            # Get the list of possible columns
            all_cnames = list_cnames(
                self.params[table_kw],
                axis=axis
            )

            # Select the column name
            if self.show_param(cname_kw):
                self.selectbox(
                    label="Column",
                    options=all_cnames,
                    kw=cname_kw,
                    copy_to=label_kw
                )

            # Input the column label
            if self.show_param(label_kw):
                self.text_input(label="Label", kw=label_kw)

            # If the column is a color column
            if col_elem.get("colorscale", False):


                # Select a pallete for the colors

                # Let the user select whether the colors are categorical
                if self.show_param(is_categorical_kw):
                    self.checkbox(
                        "Categorical Values",
                        kw=is_categorical_kw,
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

                    if self.show_param(colorscale_kw):
                        self.selectbox(
                            label="Select color scale",
                            options=colors_qualitative,
                            kw=colorscale_kw
                        )

                else:

                    if self.show_param(colorscale_kw):
                        self.selectbox(
                            label="Select color scale",
                            options=px.colors.named_colorscales(),
                            kw=colorscale_kw
                        )

    def render_dataframe_tables(
        self,
        key,
        axis: int
    ):
        tables_kw = join_kws(key, "tables")

        if self.mdata is None:

            all_tables = tree_tables(axis)

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
        if self.show_param(tables_kw):

            self.multiselect(
                label="Select table(s)",
                options=all_tables,
                kw=tables_kw
            )

        # Make a DataFrame with the selected table(s)
        selected_tables: List[str] = self.params.get(tables_kw, [])
        if len(selected_tables) == 0:
            st.write("No tables selected.")
            self.params_complete = False
            return

        return join_dataframe_tables(
            selected_tables,
            axis,
            mdata=self.mdata
        )

    def build_dataframe(self, key: str, columns: dict, axis: int):
        # If any columns did not have a table selected,
        # then the DataFrame cannot be built
        for col_kw, col_elem in columns.items():
            # If the column is being used
            if (
                col_elem.get("optional", False) is False
                or
                self.param(key, col_kw, "enabled")
            ):
                # The user may have selected 1 or more tables
                selected_tables = self.param(key, col_kw, "table")

                # If the column has no tables selected
                if (
                    selected_tables is None
                    or selected_tables == []
                ):
                    self.params_complete = False
                    return

        # Get the information for each column
        col_data = {
            col_kw: get_dataframe_column(
                mdata=self.mdata,
                axis=axis,
                **{
                    kw: self.param(key, col_kw, kw)
                    for kw in ["table", "cname"]
                }
            )
            for col_kw, col_elem in columns.items()
            if (
                self.param(key, col_kw, "table") is not None
                and
                len(self.param(key, col_kw, "table")) > 0
            )
            if (
                col_elem.get("optional", False) is False
                or
                self.param(key, col_kw, "enabled")
            )
        }

        # Build a DataFrame with no null values
        return pd.DataFrame(col_data).dropna()

    def filter_dataframe_cols(
        self,
        key: str,
        axis: int,
        df: pd.DataFrame
    ) -> pd.DataFrame:
        """
        Let the user select a subset of columns for analysis.
        """
        cols_query_kw = join_kws(key, "cols_query")
        if self.show_param(cols_query_kw):
            if self.ix == -1:
                container = st.expander("Filter Columns")
            else:
                container = st.container(border=True)
                container.write("**Filter Columns**")
        else:
            container = st.container()

        with container:
            df = self.render_query(
                cols_query_kw,
                0 if axis else 1,  # Filter the opposite axis
                1,  # Perform the filtering on the columns
                df
            )

        if self.show_param(cols_query_kw) and df is not None:
            container.write(f"Number of filtered columns: {df.shape[1]:,}")

        return df

    def filter_dataframe_rows(
        self,
        key: str,
        axis: int,
        df: pd.DataFrame
    ) -> pd.DataFrame:
        """
        Let the user select a subset of rows for analysis.
        """
        rows_query_kw = join_kws(key, "rows_query")
        if self.ix == -1:
            container = st.expander("Filter Rows")
        elif self.show_param(rows_query_kw):
            container = st.container(border=True)
            container.write("**Filter Rows**")
        else:
            container = st.container()

        with container:
            df = self.render_query(
                rows_query_kw,
                axis,
                0,  # Perform the filtering on the rows
                df
            )

        if self.show_param(rows_query_kw) and df is not None:
            container.write(f"Number of filtered rows: {df.shape[0]:,}")
        return df

    def render_query(
        self,
        key: str,
        axis: int,
        filter_axis: int,
        df: pd.DataFrame
    ) -> pd.DataFrame:

        type_kw = join_kws(key, "query", "type")

        # Set the default values
        if self.mdata is None:

            # Type of selection, either by value or by index
            if self.params.get(type_kw) is None:
                # By default, select by value
                self.update_view_param(
                    type_kw,
                    "value"
                )

            # Get the list of tables available for all modalities
            all_tables = tree_tables(axis)

            # Table selection
            table_kw = join_kws(key, "query", "table")
            if self.params.get(table_kw) is None:
                self.update_view_param(table_kw, [])
            if any([
                table not in all_tables
                for table in self.params.get(table_kw)
            ]):
                self.update_view_param(
                    table_kw,
                    [
                        table
                        for table in self.params.get(table_kw)
                        if table in all_tables
                    ]
                )

            # Get the list of possible columns
            all_cnames = list_cnames(
                self.params[table_kw],
                axis=axis
            )

            # Column name selection
            cname_kw = join_kws(key, "query", "cname")
            if self.params.get(cname_kw) not in all_cnames:

                self.update_view_param(
                    cname_kw,
                    (
                        all_cnames[0]
                        if len(all_cnames) > 0
                        else ""
                    )
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
        if self.show_param(type_kw):

            # First figure out whether we're selecting by value
            # or by selecting specific indices
            names = [
                "Filtering by Value",
                "Selecting Specific {}".format(
                    "Rows" if filter_axis == 0 else "Columns"
                )
            ]
            values = ["value", "index"]
            self.selectbox(
                label="Select by:",
                options=values,
                display_options=names,
                kw=type_kw,
                invalidate=[value_kw]
            )

            # If we're selecting by value
            if self.params[type_kw] == "value":

                # Select the table of interest
                self.multiselect(
                    label="Table",
                    options=all_tables,
                    kw=table_kw
                )

                # If no tables have been selected
                if len(self.params.get(table_kw, [])) == 0:
                    st.write("No tables selected.")
                    return df

                # Get the list of possible columns
                all_cnames = list_cnames(
                    self.params[table_kw],
                    axis=axis
                )

                # Select the column name
                self.selectbox(
                    label="Column",
                    options=all_cnames,
                    kw=cname_kw
                )
                operators = [">=", "<=", "==", "!=", ">", "<", "in", "not in"]

            else:
                operators = ["in", "not in"]

            # Input the comparison operator
            self.selectbox(
                label="Keep items which are:",
                options=operators,
                kw=expr_kw
            )

            if self.params[expr_kw] in ["in", "not in"]:

                # If we're selecting by value
                if self.params[type_kw] == "value":

                    # Let the user select specific values, selecting
                    # from the values in the table
                    value_options = get_dataframe_column(
                        None,
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
                self.multiselect(
                    label="Value",
                    options=value_options,
                    help="Select the values for filtering",
                    kw=value_kw
                )

            else:

                # Input the comparison value directly
                self.text_input(
                    label="Value",
                    kw=value_kw,
                    help="Enter a value to filter samples on."
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
            if self.show_param(type_kw):
                st.write("Provide a value to filter samples.")
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

        # Get the table
        table = join_dataframe_tables(
            query["table"],
            axis,
            mdata=self.mdata
        )
        if table is None:
            if self.params_editable:
                st.write("No data available for filtering.")
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
                    st.write("No values match the filter criteria.")
                return df

            # Get the values from the column
            vals = table[query["cname"]]

            # If there are multiple columns with the same name
            if len(vals.shape) > 1:
                # Only take the first one
                vals = vals.iloc[:, 0]

            if query["expr"] == "in":
                return df.loc[
                    vals.isin(query["value"])
                ]
            else:
                return df.loc[
                    ~vals.isin(query["value"])
                ]

        # Apply the filter, trying both string and numeric values
        filtered_table = self._filter_table(table, query)

        # If no values are returned
        if filtered_table is None:
            st.write("No samples match the filter criteria.")
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

    def transform_dataframe(
        self,
        key: str,
        df: pd.DataFrame
    ) -> pd.DataFrame:

        # Get the list of transformations
        transforms_kw = join_kws(key, "transforms")

        if self.ix == -1:
            container = st.expander("Transform Values")
        elif self.show_param(key):
            container = st.container(border=True)
            container.write("**Transform Values**")
        else:
            container = st.container()

        # Get the list of transformations
        transforms = self.param(transforms_kw)

        if self.show_param(transforms_kw):

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
            tr = get_transform(transform)
            try:
                df = tr.run(df)
            except Exception as e:
                if container is not None:
                    container.exception(e)
                else:
                    raise e
            if self.show_param(transforms_kw):
                container.write(f"Ran: {tr.name}")

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

    @st.dialog("Add Transformation")
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

    def display(self):
        pass

    def toggle_sidebar(self, kw):
        """Show a checkbox to toggle the sidebar."""
        sidebar_kw = join_kws(kw, "sidebar")
        sidebar_param_key = self.param_key(sidebar_kw)

        st.checkbox(
            "Show in sidebar",
            self.param_in_sidebar(kw),
            key=sidebar_param_key,
            help="Optionally display this input menu item in the sidebar",
            on_change=self._toggle_sidebar,
            args=(kw,)
        )
        if st.session_state[sidebar_param_key] != self.param(sidebar_kw):
            self.update_view_param(sidebar_kw, st.session_state[sidebar_param_key])

    def _toggle_sidebar(self, kw):
        sidebar_param_key = self.param_key(join_kws(kw, "sidebar"))
        value = st.session_state[sidebar_param_key]
        # If this is different than the params
        if value != self.param(kw, "sidebar"):
            # Update the view in the mdata object
            self.update_view_param(join_kws(kw, "sidebar"), value)

    def param_in_sidebar(self, kw):
        return self.param(kw, "sidebar", default=False)

    def show_param(self, *kws: List[str]) -> bool:
        """Show a param if ?editable=True and either self.params_editable or param_in_sidebar."""
        if not get_edit_views_flag():
            return False

        if self.params_editable:
            return True
        for kw in kws:
            if self.params_editable or self.param_in_sidebar(kw):
                return True
        return False

@st.dialog("Figure Parameters", width='large')
def _show_view_sdk_snippet(ix: int):
    view = get_views()[ix]
    st.code(sdk_snippet(view))


def sdk_snippet(view: dict):
    assert "type" in view.keys()
    assert "params" in view.keys()
    params = nest_params(view["params"])
    params_str = (
        json.dumps(params, indent=4)
        .replace('false', 'False')
        .replace('true', 'True')
        .replace('null', 'None')
        .replace("\n", "\n    ")
    )
    return f"""view.{view['type'].replace('-', '_')}(
    mdata,
    **{params_str}
)
"""

