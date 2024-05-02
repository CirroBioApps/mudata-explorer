# MuData Explorer
Exploratory data analysis of multimodal datasets

## About

The MuData Explorer is a data science application developed and maintained by
the Data Core at the Fred Hutch Cancer Center.
It is primarily intended to serve the needs of Fred Hutch investigators for the
analysis of research data.

> The publicly hosted instance of the MuData Explorer should not be used for
any human data or private health information.

Any bug reports or inquiries should be directed to:
[github.com/CirroBioApps/mudata-explorer/issues](https://www.github.com/CirroBioApps/mudata-explorer/issues)

## Multimodal Data

Data can be anything or take any form.
By describing data as 'multimodal' we are proscribing a form which can be used
for many different types of data, and can be used to make analysis and visualization
easier to execute.

First break up data into observations and variables.
An observation is a thing that can be measured - a sample, person, place, thing.
A variable is an aspect of that thing which can be measured - temperature, weight,
pH, density, gene expression level, microbial abundance.

Next break up data and metadata.
This distinction is pretty arbitrary, but it helps to think about variables which
are categorical in nature, or which are available _a priori_ as being metadata.
When I collect a set of samples and I'm interested in the differences between the
treatment and control group, the treatment vs. control label is metadata while the
measurements I make about those samples are data.

When collecting data, there may be different groups of measurements which go together,
or modalities.
For example, in a particular experiment I may take observations on a daily basis
from a group of mice which were assigned to treatment and control groups.
On each day I collect a blood sample and a stool sample.
Each blood sample is processed for RNAseq to measure gene expression, as well as
mass spectrometry to assess relative abundance of different metabolites.
Each stool sample is processed to measure the relative abundance of different
bacterial species.

In this scenario, three different _modalities_ of data have been collected:

- Gene expression levels
- Metabolite abundance levels
- Microbial relative abundances

It makes sense to keep track of each modality of data separately, because the
relative expression of a particular gene is influenced by the relative expression
of all of the other genes (since the measurement is inherently proportional)
but measured totally independently from the relative abundance of any bacterial
species.

The goal of this application is to make it easier for researchers to aggregate
datasets, perform basic data science analyses, configure visualizations, and
share those analysis + visualizations with colleagues.

## Basic Features

- Add Datasets
- Add Metadata
- Perform Analysis
- Configure Visualizations
- Save and Load

## Saving and Loading Datasets

When a dataset is saved to a file, that file includes:

- Metadata (annotations)
- Data (observations)
- Analysis results
- Analysis history
- Visualizations
- Unique hash of data contents

If you load a file that has been provided by a colleague, it will automatically
load any visualizations ("Views") which have been set up, as well as the data used
to generate that visualization.

## Hashing Datasets

Whenever a dataset is saved to a file, that file's name will include a "hash."
That hash is a unique string which is generated using the SHA256 library as a
function of the total content of that dataset (including the history and visualizations).
When that dataset is loaded into the app, the hash in the file name will be checked
against the file contents and a message will be displayed if those values match.
If a file is modified in transit or at rest but the filename remains constant, it will
be immediately obvious when the file is loaded.

### Known Bugs -- Hashing

When a dataset contains null values (i.e. any missing values in the data or metadata),
the hash will change when the dataset is saved to a file for the first time.
The root cause of this bug has to do with how null values are treated in Python.
To work around this bug for any datasets which contain null values, simply perform the
save/load cycle one extra time before sharing the file with a colleague.

## Developers Notes

### Shared State

The entirety of the app is driven by the shared `st.session_state['mdata']` element,
which contains a single `muon.MuData` object.
Each tab in the app is set up as a separate page which loads and operates independently.
Any changes made on any particular page are saved by modifying elements of the `MuData`
object, specifically in the `.uns` slot.

- `mudata-explorer-views`: A `list` of `dicts` with all of the data needed to render each visualization
- `mudata-explorer-settings`: A `dict` with any app-wide settings
- `mudata-explorer-history`: A `list` of each event recorded in the history for the dataset
- `mudata-explorer-provenance`: A `dict` describing the event which resulted in the creation of any data element
