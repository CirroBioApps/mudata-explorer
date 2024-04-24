# mudata-explorer
Exploratory data analysis of multimodal datasets

## Developers Notes

### Data Views

To provide a flexible interface for data exploration, each display element of
the mudata-explorer will be configured as a separate "View" class.
Each View will be driven by a set of parameters and will have methods for:

1. Processing
2. Display

### Order of Operations

1. Check if the `processed` flag has been set to `true` on the View:
    - If `true`, skip to step 3
    - If `false`, continue to step 2
2. When the user clicks the "Run" button, execute the `Processing` method
    - Set the `processed` flag to `true`
    - Capture the logs from the `Processing` method
3. Display the logs from the `Processing` method
4. Execute the `Display` method

### Data Model

The core data necessary to drive the entire visualization consists of:

1. An AnnData object
2. A list of dicts containing the params provided to each View

> Note: The list of views is explicitly ordered to ensure that the views are
executed in the correct order.
