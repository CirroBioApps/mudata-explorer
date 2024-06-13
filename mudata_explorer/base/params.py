from dataclasses import dataclass
from typing import Dict, List, Optional, Any, Union
from streamlit.delta_generator import DeltaGenerator


@dataclass
class Param:
    _key: str  # Unique string that identifies this element
    is_complete: bool
    streamlit_enabled: bool
    streamlit_container: Optional[DeltaGenerator]

    def key(self, *args):
        return ".".join([
            i
            for i in [self._key, *args]
            if i is not None
        ])


@dataclass
class StringParam(Param):
    default: Optional[str]
    optional: bool
    enabled: bool


@dataclass
class FloatParam(Param):
    default: Optional[str]
    optional: bool
    enabled: bool


@dataclass
class BooleanParam(Param):
    default: Optional[str]
    optional: bool
    enabled: bool


@dataclass
class IntegerParam(Param):
    default: Optional[str]
    optional: bool
    enabled: bool


@dataclass
class SupportingFigureParam(Param):
    default: Optional[str]
    optional: bool
    enabled: bool


@dataclass
class ColorscaleParam(Param):
    is_categorical: bool
    scale: Optional[str]


@dataclass
class ColumnParam(Param):
    tables: List[str]
    cname: str
    label: str
    optional: bool
    enabled: bool
    colorscale: Optional[ColorscaleParam]


@dataclass
class RowsQueryParam(Param):
    type: str  # 'value' or 'index'
    tables: List[str]
    cname: str
    expr: str
    value: str


@dataclass
class ColsQueryParam(Param):
    type: str  # 'value' or 'index'
    tables: List[str]
    cname: str
    expr: str
    value: str


@dataclass
class DataFrameColumnsParam(Param):
    """Build a DataFrame from a list of columns."""
    axis: int
    columns: List[ColumnParam]
    rows_query: RowsQueryParam
    cols_query: ColsQueryParam
    transforms: List[str]


@dataclass
class DataFrameTablesParam(Param):
    """Build a DataFrame from a list of tables."""
    axis: int
    tables: List[str]
    rows_query: RowsQueryParam
    cols_query: ColsQueryParam
    transforms: List[str]


class Object:
    def __init__(
        self,
        params: Dict[str, Union[Param, 'Object']],
        label: Optional[str] = None,
    ):
        self._params = params

    @property
    def is_complete(self):
        return all([
            param.is_complete
            for param in self._params.values()
        ])


def hydrate_param(
    prefix: str,
    elem: dict,
    **kwargs
):
    assert "type" in elem
    if elem["type"] == "object":
        return Object(
            params={
                key: hydrate_param(
                    (
                        f"{prefix}.{key}"
                        if prefix is not None
                        else key
                    ),
                    sub_elem,
                    **kwargs
                )
                for key, sub_elem in elem["properties"].items()
            },
            label=elem.get("label")
        )
    elif elem["type"] == "string":
        return StringParam(
            _key=prefix,
            default=elem.get("default"),
            optional=elem.get("optional", False),
            enabled=elem.get("enabled", True),
            **kwargs
        )

    elif elem["type"] == "float":
        return FloatParam(
            _key=prefix,
            default=elem.get("default"),
            optional=elem.get("optional", False),
            enabled=elem.get("enabled", True),
            **kwargs
        )

    elif elem["type"] == "boolean":
        return BooleanParam(
            _key=prefix,
            default=elem.get("default"),
            optional=elem.get("optional", False),
            enabled=elem.get("enabled", True),
            **kwargs
        )

    elif elem["type"] == "integer":
        return IntegerParam(
            _key=prefix,
            default=elem.get("default"),
            optional=elem.get("optional", False),
            enabled=elem.get("enabled", True),
            **kwargs
        )

    elif elem["type"] == "supporting_figure":
        return SupportingFigureParam(
            _key=prefix,
            default=elem.get("default"),
            optional=elem.get("optional", False),
            enabled=elem.get("enabled", True),
            **kwargs
        )

    elif elem["type"] == "colorscale":
        return ColorscaleParam(
            _key=prefix,
            is_categorical=elem.get("is_categorical", False),
            scale=elem.get(
                "scale",
                (
                    "D3"
                    if elem.get("is_categorical", False)
                    else "Viridis"
                )
            ),
            **kwargs
        )

    elif elem["type"] == "column":
        return ColumnParam(
            _key=prefix,
            tables=elem["tables"],
            cname=elem["cname"],
            label=elem["label"],
            optional=elem.get("optional", False),
            enabled=elem.get("enabled", True),
            colorscale=elem.get("colorscale"),
            **kwargs
        )

    elif elem["type"] == "rows_query":
        return RowsQueryParam(
            _key=prefix,
            type=elem["type"],
            tables=elem["tables"],
            cname=elem["cname"],
            expr=elem["expr"],
            value=elem["value"],
            **kwargs
        )

    elif elem["type"] == "cols_query":
        return ColsQueryParam(
            _key=prefix,
            type=elem["type"],
            tables=elem["tables"],
            cname=elem["cname"],
            expr=elem["expr"],
            value=elem["value"],
            **kwargs
        )

    elif elem["type"] == "dataframe_columns":
        return DataFrameColumnsParam(
            _key=prefix,
            axis=elem["axis"],
            columns=[
                hydrate_param(
                    f"{prefix}.columns.{i}",
                    column,
                    **kwargs
                )
                for i, column in enumerate(elem["columns"])
            ],
            rows_query=hydrate_param(
                f"{prefix}.rows_query",
                elem["rows_query"],
                **kwargs
            ),
            cols_query=hydrate_param(
                f"{prefix}.cols_query",
                elem["cols_query"],
                **kwargs
            ),
            transforms=elem.get("transforms", []),
            **kwargs
        )

    elif elem["type"] == "dataframe_tables":
        return DataFrameTablesParam(
            _key=prefix,
            axis=elem["axis"],
            tables=elem["tables"],
            rows_query=hydrate_param(
                f"{prefix}.rows_query",
                elem["rows_query"],
                **kwargs
            ),
            cols_query=hydrate_param(
                f"{prefix}.cols_query",
                elem["cols_query"],
                **kwargs
            ),
            transforms=elem.get("transforms", []),
            **kwargs
        )

    else:
        raise ValueError(f"Unknown type: {elem['type']}")


class Schema(Object):

    pass
