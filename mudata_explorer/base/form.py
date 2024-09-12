from mudata_explorer.helpers.join_kws import join_kws
from mudata_explorer.app.query_params import get_edit_view_flag
from mudata_explorer.app.session_state import get_show_sidebar_flag
from mudata_explorer.app.mdata import get_supp_figs, tree_tables, join_dataframe_tables
from mudata_explorer.app.mdata import get_view, set_view
from mudata_explorer.base import all_transforms, get_transform
import pandas as pd
import plotly.express as px
import streamlit as st
from typing import List, Optional, Dict, Any, Union


def _render_key(ix: int, prefix: str) -> str:
    return f"form-{ix}-{prefix}"


def _watch(ix: int, prefix: str, copy_to=None):
    """
    Watches for any changes to the value of the input element,
    and updates the value of the element accordingly.
    Can be overridden by child classes.
    """
    render_key = _render_key(ix, prefix)
    if render_key in st.session_state:
        value = st.session_state[render_key]
        _save_value(ix, prefix, value)
        _save_copy(ix, copy_to, value)


def _save_copy(ix, copy_to, value):
    """Save the value to a copy_to location if provided."""
    if copy_to is not None:
        if isinstance(copy_to, list):
            for prefix in copy_to:
                _save_value(ix, prefix, value)
        else:
            _save_value(ix, copy_to, value)


def _watch_enum(
    ix: int,
    prefix: str,
    enum: List[str],
    enumNames: Optional[List[str]],
    copy_to: Optional[Union[str, List[str]]] = None
):
    """Account for the case of enumNames."""
    render_key = _render_key(ix, prefix)
    value = st.session_state[render_key]
    if value is None:
        return
    if enumNames is not None:
        value = enum[enumNames.index(value)]
    _save_value(ix, prefix, value)
    _save_copy(ix, copy_to, value)


def _watch_enum_multi(
    ix: int,
    prefix: str,
    enum: List[str],
    enumNames: Optional[List[str]],
    copy_to: Optional[Union[str, List[str]]] = None
):
    """Account for the case of enumNames."""
    render_key = _render_key(ix, prefix)
    values = st.session_state[render_key]
    if values is None:
        values = []
    if enumNames is not None:
        values = [
            enum[enumNames.index(val)]
            for val in values
        ]
    _save_value(ix, prefix, values)
    _save_copy(ix, copy_to, values)


def _save_value(ix: int, prefix: str, value: Any):
    """
    Save the value of the current element in the state.
    """

    # Get the params saved in the state for this view
    view = get_view(ix)
    if "params" not in view:
        view["params"] = {}

    # The keyword for the params is {prefix}.value for everything except the sidebar
    kw = prefix if prefix.endswith(".sidebar") else f"{prefix}.value"

    # If the value does not match
    if view["params"].get(kw) != value:

        # Update it
        view["params"][kw] = value

        # Save the state
        set_view(ix=ix, view=view)


class _SharedFunctions:
    def _kw(self, *kws):
        """Join the provided kws to self.prefix"""
        return join_kws(self.prefix, *kws)
    
    def _render_key(self) -> str:
        """Unique key used to identify an input element for a single form element."""
        return _render_key(self.ix, self.prefix)
    
    def _nest_outermost(self, schema: dict, prefix: Optional[str], ix: int):
        # If we are at the top level of the schema
        if prefix is None:
            # And there is no type defined
            if "type" not in schema:
                # Assume that the outermost scope has been omitted
                return {
                    "type": "object",
                    "properties": schema
                }
        return schema

    def load(self, params: Union[Any, Dict[str, Any]]):
        """
        Load the basic values defined in a dict of params.
        """
        if not isinstance(params, dict):
            self.value = params
            return

        if "value" in params:
            self.value = params["value"]
        if "sidebar" in params and hasattr(self, 'sidebar'):
            self.sidebar.value = params["sidebar"]


class MuDataAppSidebarToggle(_SharedFunctions):
    """
    Element used to toggle whether an element will be shown in the sidebar.
    """

    label = "Show in Sidebar"
    value: bool

    def __init__(self, value: bool, prefix: str, ix: int):
        self.prefix = join_kws(prefix, "sidebar")
        self.ix = ix
        self.value = st.session_state.get(self._render_key(), value)

    def render(self, copy_to=None):
        """Display the checkbox and update the value on change."""
        st.checkbox(
            "Show in Sidebar",
            value=self.value,
            key=self._render_key(),
            on_change=_watch,
            args=(self.ix, self.prefix,),
            kwargs=dict(copy_to=copy_to)
        )


class MuDataAppEnabledToggle(_SharedFunctions):
    """
    Element used to toggle whether an element is enabled.
    """

    label = "Enabled"
    value: bool
    # Whether or not the toggle element is displayed in the sidebar
    sidebar: MuDataAppSidebarToggle

    def __init__(self, value: bool, prefix: str, ix: int, label=None):
        self.prefix = join_kws(prefix, "enabled")
        self.ix = ix
        self.value = st.session_state.get(self._render_key(), value)
        self.sidebar = MuDataAppSidebarToggle(
            value=False,
            prefix=join_kws(prefix, "sidebar"),
            ix=ix
        )
        if label is not None:
            self.label = label

    def render(self, copy_to=None):
        """Display the checkbox and update the value on change."""
        st.checkbox(
            self.label,
            value=self.value,
            key=self._render_key(),
            on_change=_watch,
            args=(self.ix, self.prefix,),
            kwargs=dict(copy_to=copy_to)
        )
        self.value = st.session_state[self._render_key()]


class MuDataAppDummy:
    """
    Dummy element used instead of MuDataAppEnabledToggle
    when an element is not optional,
    or instead of MuDataAppSidebarToggle when an element
    cannot be shown in the sidebar.
    """
    value = False

    def __init__(self):
        return

    def render(self, copy_to=None):
        return

    def dump(self, mdata=None) -> Dict[str, Any]:
        return {}
    
    def dehydrate(self):
        return {}
    
    def load(self, _):
        return


class MuDataAppFormElement(_SharedFunctions):
    """Attributes and functions shared by all elements of a form."""

    # ix: the index of the top-level Action (View / Process)
    ix: int
    # prefix: the path to the element (e.g. formatting.max_features.value)
    prefix: str
    # Human readable label
    label: str
    # Optional help text
    help: Optional[str]
    # Value of the element
    value: Any
    # Whether or not the form element is displayed in the sidebar
    sidebar: Union[MuDataAppSidebarToggle, MuDataAppDummy]
    # Whether the form element can be toggled to null value
    optional: bool
    # If optional, use a toggle to track whether the value is null
    enabled: Union[MuDataAppEnabledToggle, MuDataAppDummy]
    # Flag used to indicate whether sufficient input has been provided
    complete: bool

    def __init__(self, schema: dict, prefix=None, ix=-1, sidebar=True):
        """
        Set up the form element.
        The schema defines the attributes of the element.
        The prefix defines its position nested within a
        larger form.
        The ix indicates the index position of the form
        within the app as a whole.
        """

        assert isinstance(schema, dict)
        schema = self._nest_outermost(schema, prefix, ix)
        assert "type" in schema, schema

        self.prefix = prefix
        self.ix = ix
        self.label = schema.get(
            "label",
            prefix.split(".")[-1] if isinstance(prefix, str) else None
        )
        self.help = schema.get("help")
        self.value = st.session_state.get(self._render_key(), schema.get("default"))
        if sidebar:
            self.sidebar = MuDataAppSidebarToggle(schema.get("sidebar", False), prefix, ix)
        else:
            self.sidebar = MuDataAppDummy()

        self.optional = schema.get("optional", False)
        self.enabled = MuDataAppEnabledToggle(
            schema.get("enabled", True),
            prefix,
            ix,
            label=f"{self.label}: Enabled"
        )

        self._complete = True

    @property
    def complete(self):
        return self._complete

    @property
    def enabled_or_required(self) -> bool:
        """If the element is either not optional or enabled."""
        return self.enabled.value if self.optional else True

    @property
    def enabled_or_required_in_sidebar(self) -> bool:
        """
        If the element is either not optional or enabled, while also
        marked for being shown in the sidebar.
        """
        if self.sidebar.value:
            return self.enabled.value if self.optional else True
        else:
            return False
    
    @property
    def enabled_or_required_in_sidebar_recur(self) -> bool:
        """Dummy method overridden by form elements with sub elements."""
        return self.enabled_or_required_in_sidebar
    
    @property
    def show_in_sidebar(self) -> bool:
        """Whether the element should be shown in the sidebar."""
        return self.sidebar.value
    
    @property
    def show_in_sidebar_recur(self) -> bool:
        """
        Whether the element should be shown in the sidebar,
        recursing into any child elements.
        """
        return self.sidebar.value

    @property
    def _in_sidebar(self) -> bool:
        """Wrapper to determine if the form is being rendered in the sidebar."""
        return self.ix >= 0 and get_edit_view_flag() != self.ix and get_show_sidebar_flag()

    @property
    def _display_only(self) -> bool:
        """Wrapper to determine if no form inputs should be shown."""
        return self.ix >= 0 and get_edit_view_flag() != self.ix and not get_show_sidebar_flag()

    def dump(self, mdata=None) -> Dict[str, Any]:
        """
        Return a dict with all values defined by the form element.
        Note that the .value attribute is also being returned as
        using the prefix key for the element overall.
        This accomplishes the back-compatibility with versions
        that were saved by the earliest versions of the app.
        """
        vals = {
            self._kw("sidebar"): self.sidebar.value,
            self.prefix: self.value
        }
        if self.enabled_or_required:
            vals[self._kw("value")] = self.value
        if self.optional:
            vals = {
                self._kw("enabled", "value"): self.enabled.value,
                self._kw("enabled", "sidebar"): self.enabled.sidebar.value,
                **vals,
            }
        return vals
    
    def dehydrate(self) -> Dict[str, Any]:
        """
        Return a dict with only those values defined by the form
        element which are JSON serializable.
        This alternate function is needed for complex elements
        like the pandas DataFrame, which needs to be included
        in the output of dump() but should not be saved to a JSON
        string.
        Overridden by child elements.
        """
        vals = {
            self._kw("value"): self.value,
            self._kw("sidebar"): self.sidebar.value
        }
        if self.optional:
            vals = {
                self._kw("enabled", "value"): self.enabled.value,
                self._kw("enabled", "sidebar"): self.enabled.sidebar.value,
                **vals,
            }
        return vals
    
    def dehydrate_no_value(self) -> Dict[str, Any]:
        vals = {
            self._kw("sidebar"): self.sidebar.value
        }
        if self.optional:
            vals = {
                self._kw("enabled", "value"): self.enabled.value,
                self._kw("enabled", "sidebar"): self.enabled.sidebar.value,
                **vals,
            }
        return vals

    def load(self, params: Union[Any, Dict[str, Any]]):
        """
        Load all values defined in a dict of params.
        """
        # Load the 'value' and 'sidebar' elements
        super().load(params)

        if not isinstance(params, dict):
            return

        # Load the 'enabled' element
        if "enabled" in params:
            self.enabled.load(params["enabled"])

    def render(self, copy_to=None):
        """
        Present the user with a visual interface for modifying
        attributes.
        Calls the _render() method of the child class for all
        type-specific functionality
        """

        # If self.ix is -1, then this is a process and we should show all options
        if self.ix == -1:
            if self.optional:
                self.enabled.render()
                if not self.enabled.value:
                    return
            self._render(copy_to=copy_to)
        # If we are editing a single view, show all options
        elif get_edit_view_flag() == self.ix:
            # Make two columns to include the sidebar toggle
            cols = st.columns([2, 1])
            if self.optional:
                with cols[0]:
                    self.enabled.render()
                with cols[1]:
                    self.enabled.sidebar.render()
                if not self.enabled.value:
                    return
            with cols[0]:
                self._render(copy_to=copy_to)
            with cols[1]:
                self.sidebar.render()
        # If the sidebar is open
        elif get_show_sidebar_flag():
            # And this specific element has the sidebar attribute set to True
            if self.sidebar.value:
                if self.optional:
                    self.enabled.render()
                    if not self.enabled.value:
                        return
                # Then show the input element
                self._render(copy_to=copy_to)
            else:
                return
        else:
            return

    def _render(self, copy_to=None):
        """
        Just show the label.
        Can be overridden by child classes.
        """

        if self.label is not None:
            st.markdown(self.label)
        if self.help is not None:
            st.markdown(self.help)

    def parse_elem(self, elem: dict, elem_key: str):
        """Parse a single element of the form."""

        assert isinstance(elem, dict), f"Expected dict, got {type(elem)}"
        assert "type" in elem, f"Missing 'type' key in schema: {elem}"

        if elem["type"] == "object":
            return MuDataAppForm(elem, prefix=elem_key, ix=self.ix)

        elif elem["type"] == "string":
            return MuDataAppString(elem, prefix=elem_key, ix=self.ix)

        elif elem["type"] == "enum":
            return MuDataAppEnum(elem, prefix=elem_key, ix=self.ix)

        elif elem["type"] == "enum_multi":
            return MuDataAppEnumMulti(elem, prefix=elem_key, ix=self.ix)

        elif elem["type"] == "float":
            return MuDataAppFloat(elem, prefix=elem_key, ix=self.ix)
        
        elif elem["type"] == "boolean":
            return MuDataAppBoolean(elem, prefix=elem_key, ix=self.ix)
        
        elif elem["type"] == "integer":
            return MuDataAppInteger(elem, prefix=elem_key, ix=self.ix)
        
        elif elem["type"] == "supporting_figure":
            return MuDataAppSupportingFigure(elem, prefix=elem_key, ix=self.ix)
        
        elif elem["type"] == "dataframe":
            return MuDataAppDataFrame(elem, prefix=elem_key, ix=self.ix)
        
        else:
            raise ValueError(f"Unrecognized type: {elem['type']}")


class MuDataAppForm(MuDataAppFormElement):
    """
    The form element is used by the 'object' type.
    Form schemas must always be initialized with this type
    at the root.
    """

    # Elements nested within the object
    properties: Dict[str, MuDataAppFormElement]

    def __init__(self, schema: dict, prefix=None, ix=-1):

        schema = self._nest_outermost(schema, prefix, ix)
        super().__init__(
            schema=schema,
            prefix=prefix,
            ix=ix,
            sidebar=False
        )        
        assert schema["type"] == "object"

        self.properties = {
            elem_key: self.parse_elem(elem, self._kw(elem_key))
            for elem_key, elem in schema.get("properties", {}).items()
        }

    @property
    def complete(self):
        """A form is complete if all child elements are complete."""
        return all([elem.complete for elem in self.properties.values()])

    def dump(self, mdata=None) -> Dict[str, Any]:
        """
        Return a dict with all values defined by the form element,
        along with all of its nested properties.
        """
        vals = {
            kw: val
            for elem in self.properties.values()
            for kw, val in elem.dump(mdata=mdata).items()
        }
        if self.optional:
            vals = {
                self._kw("enabled", "value"): self.enabled.value,
                self._kw("enabled", "sidebar"): self.enabled.sidebar.value,
                **vals,
            }
        return vals

    def dehydrate(self) -> Dict[str, Any]:
        """
        Return a dict with all values defined by the form element,
        along with all of its nested properties.
        """
        vals = {
            kw: val
            for elem in self.properties.values()
            for kw, val in elem.dehydrate().items()
        }
        if self.optional:
            vals = {
                self._kw("enabled", "value"): self.enabled.value,
                self._kw("enabled", "sidebar"): self.enabled.sidebar.value,
                **vals,
            }
        return vals
    
    def load(self, params: Union[Any, Dict[str, Any]]):
        """
        Load all values defined in a dict of params,
        including any which map to nested properties.
        """
        super().load(params)

        if not isinstance(params, dict):
            return

        for elem_kw, elem in self.properties.items():
            if elem_kw in params:
                elem.load(params[elem_kw])

    @property
    def show_in_sidebar_recur(self) -> bool:
        """Whether the element should be shown in the sidebar."""
        return any([
            elem.show_in_sidebar_recur
            for elem in self.properties.values()
        ]) or self.show_in_sidebar

    @property
    def enabled_or_required_in_sidebar_recur(self) -> bool:
        """If any of the child elements are enabled_or_required_in_sidebar."""
        return any([
            elem.enabled_or_required_in_sidebar_recur
            for elem in self.properties.values()
        ]) or self.enabled_or_required_in_sidebar

    def render(self, copy_to=None):
        """Rendering a form will render all nested properties as well."""
        # By default, only put a border around the element if it is
        # not at the top level
        border = self.prefix is not None

        # If no element of the form is going to be rendered,
        # then we should not use a border.
        # If we are in the sidebar
        if self._in_sidebar:
            # If none of its child elements are going to be shown
            if not self.enabled_or_required_in_sidebar_recur:
                # Do not show a border
                border = False
        # Or if we are in display-only mode
        elif self._display_only:
            border = False

        # Add a border if this is not the top level
        with st.container(border=border):
            super().render(copy_to=copy_to)
            # As long as the form itself is either required or enabled
            if self.enabled_or_required:
                # Render the child elements
                for elem in self.properties.values():
                    elem.render()
    

class MuDataAppString(MuDataAppFormElement):
    value: str
    enum: Optional[List[str]]
    enumNames: Optional[List[str]]
    multiline: bool
    placeholder: str

    def __init__(self, schema: dict, prefix: str, ix: int):
        super().__init__(schema, prefix, ix)
        self.enum = schema.get("enum")
        self.enumNames = schema.get("enumNames")
        self.multiline = schema.get("multiline", False)
        self.placeholder = schema.get("placeholder", "Choose an option")
        if self.multiline:
            assert self.enum is None, "Cannot have multiline enum"
        if self.enum is None:
            assert self.enumNames is None, "Cannot have enumNames without enum"
        if self.enumNames is not None:
            assert len(self.enum) == len(self.enumNames)

    def _render(self, copy_to=None):

        kwargs = dict(help=self.help, key=self._render_key())
        if self.enum:
            # Determine the index to select
            if self.value is None:
                widget_ix = None
            else:
                if self.value not in self.enum:
                    self.value = self.enum[0]
                widget_ix = self.enum.index(self.value)
            st.selectbox(
                self.label,
                self.enum if self.enumNames is None else self.enumNames,
                index=widget_ix,
                placeholder=self.placeholder,
                on_change=_watch_enum,
                args=(self.ix, self.prefix, self.enum, self.enumNames),
                kwargs=dict(copy_to=copy_to),
                **kwargs
            )
        elif self.multiline:
            st.text_area(
                self.label,
                self.value,
                on_change=_watch,
                args=(self.ix, self.prefix,),
                kwargs=dict(copy_to=copy_to),
                **kwargs
            )
        else:
            st.text_input(
                self.label,
                self.value,
                on_change=_watch,
                args=(self.ix, self.prefix,),
                kwargs=dict(copy_to=copy_to),
                **kwargs
            )

    def update_options(self, enum, enumNames=None):
        assert isinstance(enum, list), (enum, type(enum))
        if enumNames is not None:
            assert isinstance(enumNames, list)
            assert len(enumNames) == len(enum)
            self.enumNames = enumNames
        self.enum = enum

class MuDataAppFloat(MuDataAppFormElement):
    value: float
    max_value: Optional[float]
    min_value: Optional[float]
    step: Optional[float]

    def __init__(self, schema: dict, prefix: str, ix: int):
        super().__init__(schema, prefix, ix)
        self.max_value = schema.get("max_value")
        self.min_value = schema.get("min_value")
        self.step = schema.get("step")

    def _render(self, copy_to=None):

        st.number_input(
            self.label,
            min_value=self.min_value,
            max_value=self.max_value,
            step=self.step,
            value=self.value,
            key=self._render_key(),
            on_change=_watch,
            kwargs=dict(copy_to=copy_to),
            args=(self.ix, self.prefix,)
        )

class MuDataAppBoolean(MuDataAppFormElement):
    value: bool

    def _render(self, copy_to=None):
        st.checkbox(
            self.label,
            self.value,
            key=self._render_key(),
            help=self.help,
            on_change=_watch,
            args=(self.ix, self.prefix,),
            kwargs=dict(copy_to=copy_to)
        )

class MuDataAppInteger(MuDataAppFormElement):
    value: int
    max_value: Optional[int]
    min_value: Optional[int]

    def __init__(self, schema: dict, prefix: str, ix: int):
        super().__init__(schema, prefix, ix)
        if self.value is None:
            self.value = 0
        self.value = int(self.value)
        self.max_value = schema.get("max_value")
        if self.max_value is not None:
            self.max_value = int(self.max_value)
        self.min_value = schema.get("min_value")
        if self.min_value is not None:
            self.min_value = int(self.min_value)

    def _render(self, copy_to=None):
        st.number_input(
            self.label,
            min_value=self.min_value,
            max_value=self.max_value,
            value=self.value,
            key=self._render_key(),
            help=self.help,
            on_change=_watch,
            args=(self.ix, self.prefix,),
            kwargs=dict(copy_to=copy_to)
        )


class MuDataAppEnum(MuDataAppFormElement):
    """Enum element used in multiple places to select a single option."""
    value: str
    enum: List[str]
    enumNames: Optional[List[str]]
    placeholder: str

    def __init__(self, schema: dict, prefix: str, ix: int):
        super().__init__(schema, prefix, ix)
        self.enum = schema.get("enum", [])
        self.enumNames = schema.get("enumNames")
        self.placeholder = schema.get("placeholder", "Choose an option")
        assert self.enum is not None, "Must specify enum in schema"
        if self.enumNames is not None:
            assert len(self.enum) == len(self.enumNames), "enumNames length must match enum"

    def _render(self, copy_to=None):

        # If the value is a list
        if isinstance(self.value, list):
            self.value = self.value[0]

        # Determine the index to select
        if self.value is None:
            widget_ix = None
        else:
            if self.value not in self.enum:
                self.value = self.enum[0]
            widget_ix = self.enum.index(self.value)

        st.selectbox(
            self.label,
            self.enum if self.enumNames is None else self.enumNames,
            index=widget_ix,
            placeholder=self.placeholder,
            help=self.help,
            key=self._render_key(),
            on_change=_watch_enum,
            args=(self.ix, self.prefix, self.enum, self.enumNames,),
            kwargs=dict(copy_to=copy_to)
        )

    def update_options(self, enum, enumNames=None):
        assert isinstance(enum, list), (enum, type(enum))
        if enumNames is not None:
            assert isinstance(enumNames, list)
            assert len(enumNames) == len(enum)
            self.enumNames = enumNames
        self.enum = enum

        
class MuDataAppEnumMulti(MuDataAppFormElement):
    """Enum element used in multiple places to select one or more options."""
    value: List[str]
    enum: List[str]
    enumNames: Optional[List[str]]
    placeholder: str

    def __init__(self, schema: dict, prefix: str, ix: int):
        super().__init__(schema, prefix, ix)
        self.enum = schema.get("enum", [])
        self.enumNames = schema.get("enumNames")
        self.placeholder = schema.get("placeholder", "Choose an option")
        assert self.enum is not None, "Must specify enum in schema"
        if self.enumNames is not None:
            assert len(self.enum) == len(self.enumNames), "enumNames length must match enum"

    def _render(self, copy_to=None):

        # If no options are selected
        if self.value is None:
            self.value = []

        # Only show options which are available
        self.value = [val for val in self.value if val in self.enum]

        # If enumNames are available, map the default to those
        if self.enumNames is None:
            default = self.value
        else:
            default = [
                self.enumNames[self.enum.index(val)]
                for val in self.value
            ]

        st.multiselect(
            self.label,
            options=self.enum if self.enumNames is None else self.enumNames,
            default=default,
            placeholder=self.placeholder,
            help=self.help,
            key=self._render_key(),
            on_change=_watch_enum_multi,
            args=(self.ix, self.prefix, self.enum, self.enumNames,),
            kwargs=dict(copy_to=copy_to)
        )

    def update_options(self, enum, enumNames=None):
        assert isinstance(enum, list), (enum, type(enum))
        if enumNames is not None:
            assert isinstance(enumNames, list)
            assert len(enumNames) == len(enum)
            self.enumNames = enumNames
        self.enum = enum


class MuDataAppSupportingFigure(MuDataAppEnum):
    value: dict

    def __init__(self, schema: dict, prefix: str, ix: int):

        # Get the list of all supporting figures
        all_figures = get_supp_figs()

        if len(all_figures) == 0:
            st.write("No supporting figures available.")
            self._complete = False
        
        schema["enum"] = all_figures
        super().__init__(schema, prefix, ix)


class MuDataAppDataFrame(MuDataAppFormElement):
    value: Optional[pd.DataFrame]
    _columns: Dict[str, 'MuDataAppDataFrameColumn']
    _tables: MuDataAppEnumMulti
    _transforms: MuDataAppEnumMulti
    _axis: MuDataAppString

    def __init__(self, elem: dict, prefix=None, ix=-1):
        """Set up the form dataframe element."""

        super().__init__(elem, prefix=prefix, ix=ix, sidebar=False)

        self._axis = MuDataAppString(
            dict(
                type="string",
                label="Select Orientation",
                enum=[0, 1],
                enumNames=["Observations", "Variables"],
                value=elem.get("axis", 0)
            ),
            prefix=self._kw("axis"),
            ix=self.ix
        )

        self._columns = {
            col_kw: MuDataAppDataFrameColumn(
                col_elem,
                prefix=self._kw("columns", col_kw),
                ix=self.ix,
                axis=self._axis.value
            )
            for col_kw, col_elem in elem.get("columns", {}).items()
        }

        # Element which will let the user select the tables to input
        self._tables = MuDataAppEnumMulti(
            dict(
                type="enum",
                enum=elem.get("tables", []),
                default=elem.get("tables", []),
                label="Select Table(s)"
            ),
            self._kw("tables"),
            self.ix
        )
        
        # Element which will let the user select data transformations
        self._transforms = MuDataAppEnumMulti(
            dict(
                type="enum",
                enum=[
                    transform.id
                    for transform in all_transforms().values()
                ],
                enumNames=[
                    transform.name
                    for transform in all_transforms().values()
                ],
                default=elem.get("transforms", []),
                label="Optional Data Transformations"
            ),
            self._kw("transforms"),
            self.ix
        )

        # Element which lets the user filter rows
        self._filter_rows = MuDataAppDataFrameFilterAxis(
            dict(
                label="Filter Rows",
                axis=0
            ),
            self._kw("filter_rows"),
            self.ix
        )

        # Element which lets the user filter columns
        self._filter_cols = MuDataAppDataFrameFilterAxis(
            dict(
                label="Filter Columns",
                axis=1
            ),
            self._kw("filter_cols"),
            self.ix
        )

    def _elements(self) -> List[MuDataAppFormElement]:
        """Return all of the child elements."""
        return [
            self._tables,
            self._transforms,
            self._axis,
            *[
                col_elem
                for col_elem in self._columns.values()
            ]
        ]

    @property
    def complete(self):
        """A form is complete if all child elements are complete."""
        return all([elem.complete for elem in self._elements()])

    @property
    def enabled_or_required_in_sidebar_recur(self) -> bool:
        """If any of the child elements are enabled_or_required_in_sidebar."""
        return any([
            elem.enabled_or_required_in_sidebar_recur
            for elem in self._elements()
        ]) and self.enabled_or_required_in_sidebar

    @property
    def show_in_sidebar_recur(self) -> bool:
        """Whether the element should be shown in the sidebar."""
        return any([
            elem.show_in_sidebar_recur
            for elem in self._elements()
        ]) and self.show_in_sidebar

    def render(self, copy_to=None):
        """
        Method used to render the input elements for the DataFrame.
        """
        # By default show a border
        border = True

        # However, if we're in the sidebar
        if self._in_sidebar:
            # If none of its child elements are going to be shown
            if not self.enabled_or_required_in_sidebar_recur:
                # Do not show a border
                border = False
        # Or if we are in display-only mode
        elif self._display_only:
            border = False

        with st.container(border=border):
            super().render(copy_to=copy_to)

            # If the DataFrame is disabled
            if not self.enabled_or_required:
                return

            # Determine the axis
            self._axis.render()

            # If none was selected, default to 0
            if self._axis.value is None:
                _save_value(self.ix, self._axis.prefix, 0)
                self._axis.value = 0

            # The axes used for filtering must change
            # if the user changes the orientation
            # of the overall dataframe
            self._filter_rows.axis = self._axis.value
            self._filter_cols.axis = int(not self._axis.value)

            # Get the list of tables available for this orientation
            all_tables = tree_tables(self._axis.value)

            # Populate the list of options available for selection
            self._tables.update_options(all_tables)
            for col_elem in self._columns.values():
                col_elem.axis = self._axis.value
                col_elem._table.update_options(all_tables)

            # To filter columns and rows, list the tables from the appropriate axis
            self._filter_rows.update_tables(tree_tables(self._filter_rows.axis))
            self._filter_cols.update_tables(tree_tables(self._filter_cols.axis))

            # If the user is expected to select columns
            if len(self._columns) > 0:
                # Render those options
                for col_elem in self._columns.values():
                    col_elem.render()
                    # If the column is enabled or required, but no value is provided
                    if col_elem.enabled_or_required and col_elem.value is None:
                        # Then the DataFrame is not complete
                        self._complete = False

                # Build a DataFrame with the selected columns
                self._build_dataframe_from_columns()

                # Let the user filter the rows
                self._filter_rows.render()

                # Optionally filter the rows
                if self._filter_rows.value is not None:

                    # Modify self.value in place
                    self._filter_rows_in_place()

                    if border:
                        st.write(
                            f"Rows after filtering: {self.value.shape[0]:,}"
                        )

            # Otherwise
            else:
                # Render the option to select table(s) for input
                self._tables.render()
                if self._tables.value is None or len(self._tables.value) == 0:
                    self._complete = False
                    return

                # Get the DataFrame(s) selected by the user
                self.value = join_dataframe_tables(self._tables.value, self._axis.value)

                # Let the user filter the rows and the columns
                self._filter_rows.render()
                self._filter_cols.render()

                # Optionally filter the rows
                if self._filter_rows.value is not None:
                    self.value = self.value.reindex(
                        index=[
                            val
                            for val in self._filter_rows.value
                            if val in self.value.index.values
                        ]
                    )
                    if border:
                        st.write(
                            f"Rows after filtering: {self.value.shape[0]:,}"
                        )

                # Optionally filter the columns
                if self._filter_cols.value is not None:
                    self._filter_cols_in_place()
                    if border:                    
                        st.write(
                            f"Columns after filtering: {self.value.shape[0]:,}"
                        )

            # Render the option to select data transformations
            self._transforms.render()

            # Run any specified transforms
            self._run_transforms_in_place(print_status=border)

            if self.value is not None and border:
                _rows = self.value.shape[0]
                _cols = self.value.shape[1]
                _nans = self.value.isnull().sum().sum()
                _prefix = "Data table selected"
                if _rows or _cols:
                    st.write(f"{_prefix}: {_rows:,} rows and {_cols:,} columns ({_nans:,} null values)")
                else:
                    st.write(f"{_prefix}: {_rows:,} rows and {_cols:,} columns")

    def _build_dataframe_from_columns(self):
        """Build a DataFrame as a combination of columns."""

        self.value = pd.DataFrame({
            col_kw: col_elem.value
            for col_kw, col_elem in self._columns.items()
            if col_elem.value is not None
        })

    def _filter_rows_in_place(self):
        """If filtering is enabled, modify the DataFrame in self.value"""

        self.value = self.value.reindex(
            index=[
                val
                for val in self._filter_rows.value
                if val in self.value.index.values
            ]
        )

    def _filter_cols_in_place(self):
        """If filtering is enabled, modify the DataFrame in self.value"""

        self.value = self.value.reindex(
            columns=[
                val
                for val in self._filter_cols.value
                if val in self.value.columns.values
            ]
        )

    def _run_transforms_in_place(self, print_status=False):
        if self._transforms.value is not None:
            for transform_id in self._transforms.value:
                transform = get_transform(transform_id)
                if print_status:
                    st.write(f"Running: {transform.name}")
                self.value = transform.run(self.value)


    def _build_dataframe(self, mdata=None):
        """Build the DataFrame based on the parameters in the form and save as self.value"""

        # If the user is expected to select columns
        if len(self._columns) > 0:

            # Populate self.value for each of the columns
            for col_elem in self._columns.values():
                col_elem._build_dataframe(mdata=mdata)

            # Build a DataFrame with the selected columns
            self._build_dataframe_from_columns()

            # Optionally filter the rows
            if self._filter_rows.value is not None:
                self._filter_rows_in_place()

        # Otherwise
        else:
            # Stop if no tables were given
            if self._tables.value is None or len(self._tables.value) == 0:
                return

            # Get the DataFrame(s) selected by the user
            self.value = join_dataframe_tables(self._tables.value, self._axis.value, mdata=mdata)

            # Optionally filter the rows
            if self._filter_rows.value is not None:
                self._filter_rows_in_place()

            # Optionally filter the columns
            if self._filter_cols.value is not None:
                self._filter_cols_in_place()

        # Run any specified transforms
        self._run_transforms_in_place(print_status=False)

    def dump(self, mdata=None) -> Dict[str, Any]:

        items = dict()

        # If the value attribute is not populated
        if self.value is None:

            # Populate it
            self._build_dataframe(mdata=mdata)

        if self.enabled_or_required:
            items[self._kw("dataframe")] = self.value

        items = {
            **items,
            **super().dump(mdata=mdata),
            **self._axis.dump(mdata=mdata),
            **self._transforms.dump(mdata=mdata)
        }

        # If column selection is not enabled,
        # the user will select >= 1 tables
        if len(self._columns) == 0:
            items = {
                **items,
                **self._tables.dump(mdata=mdata),
                **self._filter_cols.dump(mdata=mdata)
            }

        # Add any columns specified in the schema
        else:

            for col_elem in self._columns.values():
                items = {**items, **col_elem.dump(mdata=mdata)}

        # Settings for filtering the rows
        items = {**items, **self._filter_rows.dump(mdata=mdata)}

        return items    

    def dehydrate(self) -> Dict[str, Any]:
        """
        Omit the value attribute, which is a pandas DataFrame
        and should not be serialized to JSON.
        """
        items = {
            **super().dehydrate_no_value(),
            **self._axis.dehydrate(),
            **self._transforms.dehydrate()
        }

        # If column selection is not enabled,
        # the user will select >= 1 tables
        if len(self._columns) == 0:
            items = {
                **items,
                **self._tables.dehydrate(),
                **self._filter_cols.dehydrate()
            }

        # Add any columns specified in the schema
        else:

            for col_elem in self._columns.values():
                items = {**items, **col_elem.dehydrate()}

        # Settings for filtering the rows
        items = {**items, **self._filter_rows.dehydrate()}

        return items
    
    def load(self, params: Union[Any, Dict[str, Any]]):
        """
        Perform the inverse of dehydrate, take a set of
        serialized params and save them.
        """
        super().load(params)

        if not isinstance(params, dict):
            return

        self._axis.load(params.get("axis"))        
        self._transforms.load(params.get("transforms"))
        if len(self._columns) == 0:
            self._tables.load(
                params.get(
                    "tables",
                    params.get("table")
                )
            )
            self._filter_cols.load(params.get("filter_cols"))
        else:
            for col_kw, col_elem in self._columns.items():
                col_elem.load(
                    params.get("columns", {}).get(
                        col_kw,
                        params.get(col_kw)
                    )
                )
                # Make sure to assign the loaded axis to each column element
                col_elem.axis = self._axis.value
        self._filter_rows.load(params.get("filter_rows"))


class MuDataAppDataFrameColumn(MuDataAppFormElement):
    axis: int # Updated by the parent method
    _table: MuDataAppString
    _cname: MuDataAppString
    _label: MuDataAppString
    _scale: MuDataAppString
    _colorscale: bool
    _is_categorical: MuDataAppBoolean
    value: Optional[pd.Series]

    def __init__(self, elem: Dict, prefix=None, ix=-1, axis=0):
        self.axis = axis
        elem["type"] = "column"
        super().__init__(elem, prefix, ix)

        self._table = MuDataAppEnum(
            dict(
                type="string",
                label="Select Table",
                default=elem.get("table"),
                enum=[elem.get("table")]
            ),
            self._kw("table"),
            self.ix
        )
        self._cname = MuDataAppString(
            dict(
                type="string",
                label="Select Column",
                default=elem.get("cname")
            ),
            self._kw("cname"),
            self.ix
        )
        self._label = MuDataAppString(
            dict(
                type="string",
                label="Display Label",
                default=elem.get("label")
            ),
            self._kw("label"),
            self.ix
        )

        self._colorscale = elem.get("colorscale", False)

        self._is_categorical = MuDataAppBoolean(
            dict(
                type="boolean",
                default=elem.get("is_categorical", False),
                label="Is Categorical"
            ),
            self._kw("is_categorical"),
            self.ix
        )

        self._scale = MuDataAppString(
            dict(
                type="string",
                default=elem.get("scale"),
                label="Color Scale"
            ),
            self._kw("scale"),
            self.ix
        )

    def render(self, copy_to=None):
        """
        Method used to render the input elements for the DataFrame Column.
        """
        with st.container(border=True):
            super().render(copy_to=copy_to)

            # Select a table
            self._table.render()

            if self._table.value is None \
                or len(self._table.value) == 0 \
                or any([v is None for v in self._table.value]):
                return

            # Get the DataFrame
            df = join_dataframe_tables(self._table.value, self.axis)

            # Select the column
            self._cname.update_options(list(df.columns.values))
            self._cname.render(copy_to=self._label.prefix)

            if self._cname.value is None:
                return
            
            # Save the values of the column selected
            self.value = df[self._cname.value]
            
            if self._label.value is None:
                self._label.value = self._cname.value
            self._label.render()

            if self._colorscale:
                self._is_categorical.render()
                if self._is_categorical.value:
                    # Categorical color scales
                    scale_options = [
                        cname for cname in dir(px.colors.qualitative)
                        if (
                            not cname.startswith("_")
                            and cname != "swatches"
                        )
                    ]
                else:
                    # Continuous color scales
                    scale_options = px.colors.named_colorscales()
                self._scale.update_options(scale_options)
                self._scale.render()

    def _build_dataframe(self, mdata=None):
        """Build the DataFrame based on the parameters in the form and save as self.value"""

        if self._table.value is None or len(self._table.value) == 0:
            return

        # Get the DataFrame
        df = join_dataframe_tables(self._table.value, self.axis, mdata=mdata)

        if self._cname.value is None:
            return

        # Save the values of the column selected
        self.value = df[self._cname.value]

    def dump(self, mdata=None) -> Dict[str, Any]:

        # If the value attribute has not been populated
        if self.value is None:
            # Populate it
            self._build_dataframe(mdata=mdata)

        # If the color scale has not been populated
        if self._scale.value is None:
            self._scale.value = "D3" if self._is_categorical.value else "bluered"

        return {
            **super().dump(mdata=mdata),
            **self._table.dump(mdata=mdata),
            **self._cname.dump(mdata=mdata),
            **self._label.dump(mdata=mdata),
            **self._scale.dump(mdata=mdata),
            **{self._kw("colorscale"): self._colorscale},
            **self._is_categorical.dump(mdata=mdata)
        }

    def dehydrate(self):

        return {
            **super().dehydrate_no_value(),
            **self._table.dehydrate(),
            **self._cname.dehydrate(),
            **self._label.dehydrate(),
            **self._scale.dehydrate(),
            **{self._kw("colorscale"): self._colorscale},
            **self._is_categorical.dehydrate()
        }
    
    def load(self, params: Dict[str, Any]):
        super().load(params)
        if not isinstance(params, dict):
            return
        if "colorscale" in params:
            self._colorscale = params["colorscale"]
        self._table.load(params.get("table"))
        self._cname.load(params.get("cname"))
        self._label.load(params.get("label"))
        self._scale.load(params.get("scale"))
        self._is_categorical.load(params.get("is_categorical"))

    @property
    def complete(self):
        """A form is complete if all child elements are complete."""
        if not self.enabled_or_required:
            return True
        return all([
            elem.complete
            for elem in [
                self._table,
                self._cname,
                self._label
            ]
        ])


class MuDataAppDataFrameFilterAxis(MuDataAppFormElement):
    """Filter a single axis of a DataFrame using a flexible query input."""
    axis: int
    _type: MuDataAppString
    _tables: MuDataAppString
    _cname: MuDataAppString
    _expr = MuDataAppString
    _value_enum = MuDataAppEnumMulti
    _value_str = MuDataAppString
    value = Optional[List[str]]

    value_operators = [">=", "<=", "==", "!=", ">", "<", "in", "not in"]
    seletion_operators = ["in", "not in"]

    def __init__(self, elem: Dict, prefix=None, ix=-1):
        assert "axis" in elem
        assert elem["axis"] in [0, 1]
        self.axis = elem["axis"]
        _type = elem.get("type", "value")
        elem["type"] = "form"
        elem["enabled"] = True
        super().__init__(elem, prefix, ix)
        self.value = None

        # Type of filtering to perform
        self._type = self.parse_elem(
            dict(
                type="string",
                label="Filter by:",
                value=_type,
                enum=["value", "index"],
                enumNames=[
                    "Filtering by Value",
                    f"Selecting Specific {'Rows' if self.axis == 0 else 'Columns'}"
                ]
            ),
            self._kw("type")
        )

        # Table(s) to use for filtering (if querying by value)
        self._tables = self.parse_elem(
            dict(
                type="enum_multi",
                enum=elem.get("tables", []),
                default=elem.get("tables", []),
                label="Select Table(s) for Filtering"
            ),
            self._kw("tables")
        )

        # Column to use for filtering (if querying by value)
        self._cname = self.parse_elem(
            dict(
                type="string",
                default=elem.get("cname"),
                label=f"Filter on Values from {'Row' if self.axis else 'Column'}:"
            ),
            self._kw("cname")
        )

        # Comparison expression for filtering (if querying by value)
        self._expr = self.parse_elem(
            dict(
                type="enum",
                default=elem.get("expr"),
                label="Comparison Operation:"
            ),
            self._kw("expr")
        )

        # When the comparison value is an enum
        self._value_enum = self.parse_elem(
            dict(
                type="enum_multi",
                label="Values",
                enum=[]
            ),
            self._kw("value_enum")
        )

        # When the comparison value is a string
        self._value_str = self.parse_elem(
            dict(
                type="string",
                label="Value",
                help="Fill in the value which will be used for the comparison"
            ),
            self._kw("value_str")
        )

    def update_tables(self, all_tables: List[str]):
        """Fill in the list of tables available."""
        self._tables.update_options(all_tables)

    def render(self, copy_to=None):
        # By default show a border
        border = True

        # However, if we're in the sidebar
        if self._in_sidebar:
            # If none of its child elements are going to be shown
            if not self.enabled_or_required_in_sidebar_recur:
                # Do not show a border
                border = False
        # Or if we are in display-only mode
        elif self._display_only:
            border = False

        with st.container(border=border):
            super().render(copy_to=copy_to)

            # Prompt for the type of querying
            self._type.render()
            if self._type.value is None:
                return

            # Let the user select a table
            self._tables.render()

            # If no tables are selected, stop
            if self._tables.value is None or len(self._tables.value) == 0:
                return
            
            # Get the DataFrame selected by the user
            df = join_dataframe_tables(self._tables.value, self.axis)
            
            # If the value-based querying is selected
            if self._type.value == "value":

                # Update the column selector with all possible options
                self._cname.update_options(list(df.columns.values))

                # Let the user select a column
                self._cname.render()

                if self._cname.value is None:
                    return

                # Fill in the list of possible options for the expression
                self._expr.update_options(self.value_operators)

            # If the index selection based querying is selected
            else:
                # Only let the user select 'in' or 'not in'
                self._expr.update_options(self.seletion_operators)

            # Let the user select what type of expression to use
            self._expr.render()

            if self._expr.value is None:
                return

            # If the user opted for "in" or "not in"
            if self._expr.value in ["in", "not in"]:
                # If the value-based querying is selected
                if self._type.value == "value":
                    # Update with the list of all possible values
                    # in the selected column
                    value_options = list(df[self._cname.value].unique())

                else:
                    # Update with the list of all index values
                    # for the selected table
                    value_options = list(df.index.values)

                # Update the selector with the options as appropriate
                self._value_enum.update_options(value_options)

                # Let the user select from those options
                self._value_enum.render()

                if self._value_enum.value is None or len(self._value_enum.value) == 0:
                    return

                # Save the index of those items which match the filter

                # If the value-based querying is selected
                if self._type.value == "value":
                    if self._expr.value == "in":
                        self.value = df.index.values[df[self._cname.value].apply(lambda v: v in self._value_enum.value)]
                    else:
                        self.value = df.index.values[df[self._cname.value].apply(lambda v: v not in self._value_enum.value)]

                # If direct selection of indices was performed
                else:
                    if self._expr.value == "in":
                        self.value = self._value_enum.value
                    else:
                        self.value = [
                            val for val in value_options
                            if val not in self._value_enum.value
                        ]

            # If the user opted for a boolean operator
            else:
                # Render the element for arbitrary comparison
                self._value_str.render()

                if self._value_str.value is None or len(self._value_str.value) == 0:
                    return

                # Save the index of those items which match the filter
                try:
                    self.value = df.query(f"`{self._cname.value}` {self._expr.value} {self._value_str.value}").index.values
                except:
                    self.value = df.query(f"`{self._cname.value}` {self._expr.value} '{self._value_str.value}'").index.values

            if self.value is not None and border:
                # Tell the user how many element passed the filter
                st.write(f"Elements passing filter: {len(self.value):,}")

    def dump(self, mdata=None) -> Dict[str, Any]:
        return {
            **super().dump(mdata=mdata),
            **self._type.dump(mdata=mdata),
            **self._tables.dump(mdata=mdata),
            **self._cname.dump(mdata=mdata),
            **self._expr.dump(mdata=mdata),
            **self._value_enum.dump(mdata=mdata),
            **self._value_str.dump(mdata=mdata)
        }
        
    def dehydrate(self):
        return {
            **super().dehydrate_no_value(),
            **self._type.dehydrate(),
            **self._tables.dehydrate(),
            **self._cname.dehydrate(),
            **self._expr.dehydrate(),
            **self._value_enum.dehydrate(),
            **self._value_str.dehydrate()
        }
    
    def load(self, params: Dict[str, Any]):
        super().load(params)
        if not isinstance(params, dict):
            return
        self._type.load(params.get("type"))
        self._tables.load(params.get("tables"))
        self._cname.load(params.get("cname"))
        self._expr.load(params.get("expr"))
        self._value_enum.load(params.get("value_enum"))
        self._value_str.load(params.get("value_str"))
