# mudata-explorer
Exploratory data analysis of multimodal datasets

## Developers Notes

### Data Views

To provide a flexible interface for data exploration, each display element of
the mudata-explorer will be configured as a separate "View" class.
Each View will be driven by a set of parameters and will have methods for:

1. Display
2. Inputs

### Order of Operations

1. For each View, the app will first execute the display method, and then
2. Execute the inputs method.
3. When the inputs are updated, the new params will be saved to the state and
    the whole app will refresh.

### Data Model

The core data necessary to drive the entire visualization consists of:

1. A MuData object
2. A list of dicts containing the params provided to each View

> Note: The list of views is explicitly ordered to ensure that the views are
executed in the correct order.
