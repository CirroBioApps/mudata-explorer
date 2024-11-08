# Notes: Serialization

## Alternatives

- HDF5:
    - Pro:
        - Easily portable as a single file
        - Losslessly captures the entire MuData object
    - Con:
        - Not human readable
        - Not easily parsed for web visualization
- Views JSON:
    - Pro:
        - Optimized for web visualization
    - Con:
        - Does not save the complete dataset


### Views JSON

The architecture for saving a set of MuData views as JSON uses:

- `views.layout.json` to define the visual layout of all view (figure) elements
- `views/{uid}.json` to define each individual figure

#### `views.layout.json`

A `list` of objects, each of which can be used to arrange a list
of one or more individual figures.

The `style` of each object is compatible with flexbox syntax
for arranging a collection of flex items.

> Note: It could be possible to use another layout system, in which
case the behavior switch could be keyed on the `display` attribute.

The `items` of each object will each be rendered within the container
using the layout style defined by the container.

The `style` of each `item` will be applied to the flex item configuration.

The `uid` of each `item` will be used to source the content of each figure
from the expected filepath (`views/{uid}.json`).

**Example**

```json
[
    {
        "style": {
            "display": "flex"
        },
        "items": [
            {
                "style": {},
                "uid": "abcdefg"
            }
        ]
    }
]
```

#### `views/{uid}.json`

A `dict` with information on how a specific display
element should be rendered.

The `type` attribute is used to indicate how the `contents` should be parsed.

**Example**

```json
{
    "type": "plotly",
    "content": {"data": [], "layout": []}
}
```
