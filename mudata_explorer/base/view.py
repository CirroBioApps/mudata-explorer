from typing import Union
import pandas as pd
import plotly.express as px
import streamlit as st
from typing import List
from mudata_explorer import app
from mudata_explorer.base.base import MuDataAppHelpers
from streamlit.delta_generator import DeltaGenerator


class View(MuDataAppHelpers):

    ix: int
    type: str
    editable: bool
    name: str
    desc: str
    categories: List[str]
    params: dict = {}
    defaults: dict = {}
    view_container: DeltaGenerator
    inputs_container: DeltaGenerator
    schema: dict = {}

    def __init__(
        self,
        ix: int,
        type: str,
        editable: bool,
        name: str,
        desc: str,
        params: dict
    ):
        self.ix = ix
        self.type = type
        self.editable = editable
        self.name = name
        self.desc = desc
        self.params = {
            kw: params.get(kw, val)
            for kw, val in self.get_schema_defaults(self.schema)
        }

    def get_schema_defaults(self, schema: dict, prefix=None):
        """Yield the default values for each item in the schema."""
        for key, elem in schema.items():
            assert isinstance(elem, dict), f"Expected dict, got {type(elem)}"
            assert "type" in elem, f"Missing 'type' key in schema: {elem}"

            if elem["type"] == "object":
                yield from self.get_schema_defaults(
                    elem["properties"],
                    join_prefix_key(prefix, key)
                )

            elif elem["type"] in ["string", "number", "float", "boolean"]:
                yield (
                    join_prefix_key(prefix, key),
                    elem.get("default", None)
                )

            elif elem["type"] == "dataframe":
                for kw in ["modality", "query"]:
                    yield (
                        make_kw(prefix, key, kw),
                        elem.get(kw, None)
                    )
                if elem.get("select_columns", False):
                    yield (
                        make_kw(prefix, key, "selected_columns"),
                        []
                    )
                for col_kw, col_elem in elem.get("columns", {}).items():
                    yield (
                        make_kw(prefix, key, col_kw),
                        col_elem.get("default", None)
                    )
                    if col_elem.get("optional", False):
                        yield (
                            make_kw(prefix, key, f"{col_kw}.enabled"),
                            True
                        )
                    if col_elem.get("continuous_scale", False):
                        yield (
                            make_kw(prefix, key, f"{col_kw}.continuous_scale"),
                            "Viridis"
                        )
                    if col_elem.get("discrete_sequence", False):
                        yield (
                            make_kw(prefix, key, f"{col_kw}.discrete_sequence"),
                            "Plotly"
                        )

    def attach(self, container: DeltaGenerator):

        # Set up a container for the view
        self.view_container = container.container()

        if self.editable:
            # Set up an expander element for the parameters
            self.inputs_container = container.expander("Edit Position")

        # Now make the display, catching any errors
        try:
            self.display(self.view_container)
        except Exception as e:
            self.view_container.write(f"Error: {e}")

    @classmethod
    def template(cls):
        return dict(
            type=cls.type,
            name=cls.name,
            desc=cls.desc,
            params=cls.params
        )

    def display(self):
        """Primary method which is executed to render the view."""
        pass

    def param_key(self, kw):
        return f"view-{self.ix}-{kw}"

    def input_value_kwargs(self, kw):
        """Each input value element will be populated with default kwargs."""
        return dict(
            value=self.params[kw],
            key=self.param_key(kw),
            on_change=self.input_value_change,
            args=(kw,)
        )

    def input_selectbox_kwargs(self, kw, options: list):
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
            args=(kw,)
        )

    def input_multiselect_kwargs(self, kw):
        """Populate the multiselect element with default kwargs."""

        default = self.params[kw]

        return dict(
            key=self.param_key(kw),
            default=default,
            on_change=self.input_value_change,
            args=(kw,)
        )

    def input_value_change(self, kw):
        # Get the value provided by the user
        value = st.session_state[self.param_key(kw)]

        # If this is different than the params
        if value != self.params[kw]:

            # Update the view in the mdata object
            self.update_view_param(kw, value)

    def update_view_param(self, kw, value):
        # Get the MuData object
        mdata = app.get_mdata()

        # Modify the value of this param for this view
        mdata.uns["mudata-explorer-views"][self.ix]["params"][kw] = value

        # Save the MuData object
        app.set_mdata(mdata)

        # Also update the params object
        self.params[kw] = value

    def get_data(self, container: DeltaGenerator):

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

            prefix_key = join_prefix_key(prefix, key)

            if "label" in elem and settings["editable"]:
                container.write(f"#### {elem['label']}")

            if elem["type"] == "dataframe":

                self.render_dataframe(prefix, key, elem, container)

            elif elem["type"] == "object":

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
                            **self.input_selectbox_kwargs(prefix_key, elem["enum"])
                        )

                    else:
                        self.params[prefix_key] = container.text_input(
                            key if elem.get("label") is None else elem["label"],
                            **self.input_value_kwargs(prefix_key)
                        )

            elif elem["type"] in ["number", "float"]:

                if settings["editable"]:
                    self.params[prefix_key] = container.number_input(
                        key if elem.get("label") is None else elem["label"],
                        **self.input_value_kwargs(prefix_key)
                    )

            elif elem["type"] == "boolean":

                if settings["editable"]:
                    self.params[prefix_key] = container.checkbox(
                        key if elem.get("label") is None else elem["label"],
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

        # Set up a keyword for the modality of this dataframe
        modality_kw = make_kw(prefix, key, "modality")

        # Select the modality to use
        if settings["editable"]:
            modality_options = list(mdata.mod.keys())
            container.selectbox(
                "Select modality",
                modality_options,
                **self.input_selectbox_kwargs(
                    modality_kw,
                    modality_options
                )
            )

        # Make a single DataFrame with all of the data for this modality
        df = app.make_modality_df(mdata, self.params[modality_kw])

        # Let the user optionally filter samples
        if settings["editable"]:
            if elem.get("query") is not None:
                query_kw = make_kw(prefix, key, "query")
                container.text_input(
                    "Filter samples (optional)",
                    help="Enter a query to filter samples (using metadata or data).",
                    **self.input_value_kwargs(query_kw)
                )
                if (
                    self.params[query_kw] is not None
                    and
                    len(self.params[query_kw]) > 0
                ):
                    df = df.query(self.params[query_kw])

                    if settings["editable"]:
                        container.write(
                            f"Filtered data: {df.shape[0]:,} rows x {df.shape[1]:,} columns."
                        )

        # If there is no data, return
        if df.shape[0] == 0 or df.shape[1] == 0:
            container.write(f"No data available for {self.params['modality']}.")
            return

        if settings["editable"]:
            container.write("#### Columns")

        # Default column selection options
        cols = list(df.columns.values)

        # If the user wants to select columns
        if elem.get("select_columns", False):

            selected_cols_kw = make_kw(prefix, key, "selected_columns")

            if settings["editable"]:

                # Let the user select the columns to use
                if container.checkbox("Use all columns", value=True):
                    self.params[selected_cols_kw] = cols
                    self.update_view_param(selected_cols_kw, cols)
                else:
                    self.params[selected_cols_kw] = container.multiselect(
                        "Select columns",
                        cols,
                        **self.input_multiselect_kwargs(selected_cols_kw)
                    )

            if len(self.params[selected_cols_kw]) == 0:
                container.write("No columns selected.")
                return

            # Subset the table to just the selected columns
            cols = self.params[selected_cols_kw]
            df = df.reindex(columns=cols)

        # Iterate over the columns in the schema
        for col_kw, col_elem in elem.get("columns", {}).items():

            self.render_dataframe_column(
                prefix,
                key,
                col_kw,
                col_elem,
                cols,
                container
            )

        if len(elem.get("columns", {})) > 0:

            cols = [
                self.params[
                    make_kw(prefix, key, col_kw)
                ]
                for col_kw in elem.get("columns", {})
            ]

            cols = list(set(cols) - set([None]))

            df = df.reindex(columns=cols).dropna()

        if settings["editable"] and df.shape[0] == 0:
            container.write("No data available.")
            return

        self.params[make_kw(prefix, key, "dataframe")] = df

    def render_dataframe_column(
        self,
        prefix: Union[None, str],
        key: str,
        col_kw: str,
        col_elem: dict,
        cols: list,
        container: DeltaGenerator
    ):

        # Get the global settings
        settings = app.get_settings()

        # Get the key in the params which is used for this column
        col_path = make_kw(prefix, key, col_kw)

        # Set up a default value by picking from the list of columns
        if self.params.get(col_path) is None:
            self.update_view_param(
                col_path,
                cols.pop(0) if len(cols) > 1 else cols[0]
            )

        if settings["editable"]:
            container.selectbox(
                col_kw if "label" not in col_elem else col_elem["label"],
                cols,
                **self.input_selectbox_kwargs(col_path, cols)
            )

            # Optional column selection
            if col_elem.get("optional", False):
                enabled_kw = make_kw(prefix, key, f"{col_kw}.enabled")
                container.checkbox(
                    "Enabled",
                    **self.input_value_kwargs(enabled_kw)
                )
                if not self.params.get(enabled_kw, False):
                    self.update_view_param(col_path, None)

            # Color options
            if col_elem.get("continuous_scale", False):

                container.selectbox(
                    "Select color scale",
                    px.colors.named_colorscales(),
                    **self.input_selectbox_kwargs(
                        make_kw(prefix, key, f"{col_kw}.continuous_scale"),
                        px.colors.named_colorscales()
                    )
                )

            if col_elem.get("discrete_sequence", False):

                container.selectbox(
                    "Select color scale",
                    px.colors.named_colorscales(),
                    **self.input_selectbox_kwargs(
                        make_kw(prefix, key, f"{col_kw}.discrete_sequence"),
                        px.colors.named_colorscales()
                    )
                )


def make_kw(prefix, key, kw):
    return f"{join_prefix_key(prefix, key)}.{kw}"


def join_prefix_key(prefix, key):
    return (
        key if prefix is None else f"{prefix}.{key}"
    )
