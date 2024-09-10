from mudata_explorer.sdk import view, io, process
import pandas as pd
import mudata as mu


def prevalent_orgs_query():
    return {
        "query": {
            "cname": "prop_positive",
            "type": "value",
            "table": [
                "abund.varm.summary_stats"
            ],
            "expr": ">=",
            "value": "0.1"
        }
    }


def parse_org_label(org_label):
    """Parse the organism label into a dictionary."""

    parts = org_label.split("|")

    parts = {
        dict(
            k="kingdom",
            p="phylum",
            c="class",
            o="order",
            f="family",
            g="genus",
            s="species"
        )[x[0]]: x[3:].replace("_", " ") for x in parts
    }
    parts['index'] = parts['species']
    return parts


def read_mdata() -> mu.MuData:
    """Read in the BackhedF_2015 dataset and create a MuData object."""

    counts = pd.read_csv(
        "BackhedF_2015.relative_abundance_data.csv",
        index_col=0
    )

    # Format the taxonomic labels using the column names
    org_info = pd.DataFrame([
        parse_org_label(x)
        for x in counts.columns
    ]).set_index('index')

    # Rename the microbes to just use the species
    counts = counts.rename(
        columns=lambda x: x.split("|")[-1][3:].replace("_", " ")
    )
    meta = pd.read_csv(
        "BackhedF_2015.relative_abundance_metadata.csv",
        index_col=0
    )

    for cname in [
        "age_category",
        "pregnant",
        "family_role",
        "born_method",
        "feeding_practice"
    ]:
        meta = meta.assign(**{cname: meta[cname].fillna("n/a")})

    # Code the age metadata into broad age groups
    meta = meta.assign(
        age_group=[
            {
                "B": "A. Birth",
                "4M": "B. 4 Months",
                "12M": "C. 12 Months",
                "M": "D. Adult"
            }[sample_name.rsplit("_", 1)[1]]
            for sample_name in meta.index
        ]
    )

    # Use the abundance scaled as proportions
    abund = counts.apply(lambda x: x / x.sum(), axis=1)

    mdata = io.build_mdata(
        dict(abund=abund),
        obs=meta
    )

    # Add the org metadata
    mdata.mod['abund'].var = org_info

    return mdata


def show_header(mdata: mu.MuData):
    """Display the header of the data."""

    view.markdown(
        mdata,
        text="""# BackhedF_2015 (curatedMetagenomicData)
        
> Bäckhed F, Roswall J, Peng Y, Feng Q, Jia H, Kovatcheva-Datchary P, Li Y, Xia Y, Xie H, Zhong H, Khan MT, Zhang J, Li J, Xiao L, Al-Aama J, Zhang D, Lee YS, Kotowska D, Colding C, Tremaroli V, Yin Y, Bergman S, Xu X, Madsen L, Kristiansen K, Dahlgren J, Wang J. Dynamics and Stabilization of the Human Gut Microbiome during the First Year of Life. Cell Host Microbe. 2015 May 13;17(5):690-703. doi: 10.1016/j.chom.2015.04.004. Erratum in: Cell Host Microbe. 2015 Jun 10;17(6):852. Jun, Wang [corrected to Wang, Jun]. Erratum in: Cell Host Microbe. 2015 Jun 10;17(6):852. doi: 10.1016/j.chom.2015.05.012. PMID: 25974306.

> Abstract: 
> The gut microbiota is central to human health, but its establishment in early life has not been quantitatively and functionally examined. Applying metagenomic analysis on fecal samples from a large cohort of Swedish infants and their mothers, we characterized the gut microbiome during the first year of life and assessed the impact of mode of delivery and feeding on its establishment. In contrast to vaginally delivered infants, the gut microbiota of infants delivered by C-section showed significantly less resemblance to their mothers. Nutrition had a major impact on early microbiota composition and function, with cessation of breast-feeding, rather than introduction of solid food, being required for maturation into an adult-like microbiota. Microbiota composition and ecological network had distinctive features at each sampled stage, in accordance with functional maturation of the microbiome. Our findings establish a framework for understanding the interplay between the gut microbiome and the human body in early life.

### Sample Groups

The number of samples is shown in the table below, broken out
by the different annotations that were provided for age group,
family role, born method, and feeding practice.
        """ # noqa
    )


def show_sample_groups(mdata: mu.MuData):

    view.group_freq(
        mdata,
        **{
            "data": {
                "table": {
                    "tables": [
                        "Observation Metadata"
                    ],
                    "cols_query": {
                        "query": {
                            "cname": "",
                            "type": "index",
                            "expr": "in",
                            "value": [
                                "age_group",
                                "family_role",
                                "born_method",
                                "feeding_practice"
                            ]
                        }
                    }
                }
            },
            "options": {
                "sort": "Groupings"
            }
        }
    )

    view.markdown(
        mdata,
        text="""
The largest amount of data is available for the comparison of
samples taken at birth, 4 months, and 12 months of age.
Because a large fraction of the participants were breastfed
exclusively at birth and 4 months, there is not a large
amount of analytical power available for comparing the
effects of breastfeeding versus formula feeding within each
timepoint.
"""
    )


def show_org_summary(mdata):

    # First summarize the relative abundance of
    # the microbes across all samples
    process.summary_stats(
        mdata,
        outputs_dest_key="summary_stats",
        table_data_axis=1,
        table_data_tables=["abund.data"]
    )

    view.markdown(mdata, text="### Detected Organisms")

    # Plot the relative abundance of the microbes
    view.plotly_scatter(
        mdata,
        **{
            "data": {
                "x": {
                    "table": [
                        "abund.varm.summary_stats"
                    ],
                    "cname": "mean",
                    "label": "Avg. Abundance (Proportion)"
                },
                "y": {
                    "table": [
                        "abund.varm.summary_stats"
                    ],
                    "cname": "prop_positive",
                    "label": "Proportion of Samples Detected"
                },
                "size": {
                    "enabled": False
                },
                "color": {
                    "table": [
                        "abund.metadata"
                    ],
                    "cname": "phylum",
                    "label": "Phylum",
                    "scale": "D3",
                    "is_categorical": True
                },
                "axis": 1,
                "rows_query": prevalent_orgs_query()
            },
            "scale_options": {
                "log_x": True,
                "log_y": True
            }
        }
    )

    view.markdown(
        mdata,
        text="""The plot above summarizes the average abundance and prevalence
of the organisms that were detected in at least 10% of samples
for this dataset.
Prevalence (on the vertical axis) is defined as the proportion
of samples in which the organism was detected.
Abundance (on the horizontal axis) is defined as the proportion of the
microbial community which the organism makes up (as predicted from the
sequencing data)."""
    )

    # Show a boxplot with the most abundant organisms
    view.plotly_box_multiple(
        mdata,
        **{
            "table": {
                "category": {
                    "category": {
                        "table": [
                            "Observation Metadata"
                        ],
                        "cname": "age_group",
                        "label": "age_group"
                    },
                    "enabled": True,
                    "rows_query": {
                        "query": {
                            "cname": "age_group",
                            "type": "value",
                            "table": [
                                "Observation Metadata"
                            ],
                            "expr": "in",
                            "value": [
                                "A. Birth",
                                "D. Adult",
                                "B. 4 Months",
                                "C. 12 Months"
                            ]
                        }
                    }
                },
                "data": {
                    "tables": [
                        "abund.data"
                    ],
                    "cols_query": {
                        "query": {
                            "cname": "prop_positive",
                            "type": "value",
                            "table": [
                                "abund.varm.summary_stats"
                            ],
                            "expr": ">=",
                            "value": "0.4"
                        }
                    },
                    "rows_query": {
                        "query": {
                            "cname": "",
                            "type": "value"
                        }
                    }
                }
            },
            "variable_options": {
                "axis": "X-Axis",
                "log_values": True,
                "sort_by": "Median"
            },
            "display_options": {
                "height": 700,
                "title": "Abundant Organisms",
                "var_label": "",
                "val_label": "Relative Abundance",
                "ncols": 4
            },
            "category_options": {
                "axis": "Facet",
                "sort_by": "Name"
            }
        }
    )
    view.markdown(
        mdata,
        text="""The figure above shows the distribution of the most abundant
organisms in the dataset, with each panel showing a different age group.
Many of these abundant organisms show a striking difference in abundances
between the different age groups.
"""
    )


def workflow(
    mdata: mu.MuData,
    metadata_category: str,
    metadata_label: str,
    n_features=20,
    leiden_resolution=0.5
):
    """
    Train a classifier to predict the metadata category,
    use the weights to identify the most important features,
    run UMAP using those features, and plot the UMAP.
    Finally, cluster the samples using those features
    and show the enrichment of different clusters for the
    different age groups
    """

    # Train a classifier to predict the metadata category
    sgd_key = f"sgd_classifier_{metadata_category}"
    process.sgd_classifier(
        mdata,
        **{
            "model_params": {
                "alpha": 0.5,
                "fit_intercept": True,
                "l1_ratio": 0.15,
                "loss": "log_loss",
                "penalty": "l2",
                "random_state": 0,
                "tol": "0.001",
                "train_prop": 0.5
            },
            "outputs": {
                "dest_key": sgd_key
            },
            "table": {
                "data": {
                    "axis": 0,
                    "cols_query": prevalent_orgs_query(),
                    "tables": [
                        "abund.data"
                    ],
                    "transforms": [
                        "zscores_cols"
                    ]
                },
                "predictor": {
                    "axis": 0,
                    "predictor": {
                        "cname": metadata_category,
                        "label": metadata_category,
                        "table": [
                            "Observation Metadata"
                        ]
                    },
                    "transforms": []
                }
            }
        }
    )

    view.markdown(
        mdata,
        text=f"### Machine Learning Classifier - {metadata_label}"
    )

    view.supporting_figure(
        mdata,
        **{
            "input": "{\"attr\": \"" + sgd_key + "\", \"axis\": 0, \"modality\": \"abund\", \"slot\": \"obsm\", \"subattr\": None}:0"
        }
    )

    view.markdown(
        mdata,
        text=f"""
The figure above shows the performance of the classifier trained
to predict the {metadata_label} category.
The histogram shows the performance of the same model after
randomly permuting the labels, which provides a baseline for
evaluating the performance.
"""
    )

    view.plotly_contingency_table(
        mdata,
        **{
            "data": {
                "x": {
                    "table": [
                        f"abund.obsm.{sgd_key}"
                    ],
                    "cname": "actual",
                    "label": "actual"
                },
                "y": {
                    "table": [
                        f"abund.obsm.{sgd_key}"
                    ],
                    "cname": "predicted",
                    "label": "predicted"
                },
                "rows_query": {
                    "query": {
                        "cname": "group",
                        "type": "value",
                        "table": [
                            f"abund.obsm.{sgd_key}"
                        ],
                        "expr": "in",
                        "value": [
                            "testing"
                        ]
                    }
                }
            },
            "formatting": {
                "values": "Number of Items"
            }
        }
    )

    view.markdown(
        mdata,
        text=f"""The figure above shows the confusion matrix for 
the classifier, which compares the predicted values to the
actual values for those samples which were held out from the
training set.
This provides a detailed view of the performance of the
classifier in terms of which categories were more easily
distinguishable.

### Dimensionality Reduction and Clustering

Using the top {n_features:,} features from the classifier,
the samples were projected into a 2D space using UMAP.
In parallel, the samples were clustered using the Leiden
algorithm, which groups samples based on their similarity
in the high-dimensional space defined by the classifier's
features.
"""
    )

    # Format the query which returns a table with
    # only the most important features
    table = {
        "data": {
            "axis": 0,
            "cols_query": {
                "query": {
                    "cname": "rank",
                    "expr": "<=",
                    "table": [
                        f"abund.varm.{sgd_key}"
                    ],
                    "type": "value",
                    "value": str(n_features)
                }
            },
            "tables": [
                "abund.data"
            ]
        }
    }

    # Run UMAP using the features from the classifier
    umap_key = f"umap_{metadata_category}"
    process.umap(
        mdata,
        **{
            "outputs": {
                "dest_key": umap_key
            },
            "table": table,
            "umap_params": {
                "metric": "correlation",
                "min_dist": 0.1,
                "n_components": 2,
                "n_neighbors": 15
            }
        }
    )

    # Run leiden clustering using the same features
    leiden_key = f"leiden_{metadata_category}"
    process.leiden(
        mdata,
        **{
            "clustering": {
                "metric": "correlation",
                "n_neighbors": 15,
                "resolution": leiden_resolution
            },
            "outputs": {
                "dest_key": leiden_key
            },
            "table": table
        }
    )

    # Plot the UMAP twice, once with the metadata grouping
    # and once with the clusters
    for color in [
        {
            "cname": metadata_category,
            "is_categorical": True,
            "label": metadata_label,
            "scale": "D3",
            "table": [
                "Observation Metadata"
            ]
        },
        {
            "cname": leiden_key,
            "is_categorical": True,
            "label": f"Cluster - {metadata_label} Weighted Features",
            "scale": "D3",
            "table": [
                "Observation Metadata"
            ]
        }
    ]:
        view.plotly_scatter(
            mdata,
            **{
                "data": {
                    "color": color,
                    "size": {
                        "cname": "",
                        "enabled": False,
                        "label": "",
                        "table": []
                    },
                    "x": {
                        "cname": "UMAP 1",
                        "label": "UMAP 1",
                        "table": [
                            f"abund.obsm.{umap_key}"
                        ]
                    },
                    "y": {
                        "cname": "UMAP 2",
                        "label": "UMAP 2",
                        "table": [
                            f"abund.obsm.{umap_key}"
                        ]
                    }
                }
            }
        )

    view.markdown(
        mdata,
        text="""
The figure above shows the UMAP plot of the samples, with the
samples colored by (1) the metadata category and (2) the
clusters identified by the Leiden algorithm.
"""
    )

    # Show which clusters are enriched for different age groups
    view.plotly_contingency_table(
        mdata,
        **{
            "data": {
                "x": {
                    "table": [
                        "Observation Metadata"
                    ],
                    "cname": metadata_category,
                    "label": metadata_label
                },
                "y": {
                    "table": [
                        "Observation Metadata"
                    ],
                    "cname": leiden_key,
                    "label": f"Cluster - {metadata_label} Weighted Features"
                }
            }
        }
    )

    view.markdown(
        mdata,
        text=f"""The figure above shows the contingency table 
between the clusters assigned using the top features from the
classifier and the {metadata_label} category.
This shows which groups of samples are enriched for different
categories of {metadata_label}.

### Organisms Across Microbiome Clusters
"""
    )

    view.plotly_category_summarize_values(
        mdata,
        **{
            "table": {
                "category": {
                    "category": {
                        "table": [
                            "Observation Metadata"
                        ],
                        "cname": leiden_key,
                        "label": f"Cluster - {metadata_label} Weighted Features"
                    }
                },
                "data": {
                    "tables": [
                        "abund.data"
                    ],
                    "cols_query": {
                        "query": {
                            "cname": "rank",
                            "expr": "<=",
                            "table": [
                                f"abund.varm.{sgd_key}"
                            ],
                            "type": "value",
                            "value": str(n_features)
                        }
                    },
                    "transforms": [
                        "zscores_cols"
                    ]
                }
            },
            "formatting": {
                "colorscale": "blues",
                "sort_by": "Values",
                "color": "None"
            }
        }
    )

    view.markdown(
        mdata,
        text="""The figure above shows the average abundance of
each organism in each cluster.
To account for different scales of abundances across organisms,
the values are shown as z-scores.

### Summary

This analysis provides an example of how a machine learning approach
can be used to identify the organisms which most clearly differentiate
different groups of microbiome samples. One of the key advantages of
this approach is that it can be used to identify a collection of
different community types which may not be easily identified using
univariate analysis of individual organisms between sample groups.
"""
    )


def run():
    # Read in the data
    mdata = read_mdata()

    # Display the header
    show_header(mdata)

    # Show the number of samples from different groups
    show_sample_groups(mdata)

    # Show all of the detected organisms
    show_org_summary(mdata)

    # Train a classifier to predict the age group,
    # use the weights to identify the most important features,
    # run UMAP using those features, and plot the UMAP.
    # Finally, cluster the samples using those features
    # and show the enrichment of different clusters for the
    # different age groups
    workflow(
        mdata,
        "age_group",
        "Age Group",
        n_features=20,
        leiden_resolution=0.25
    )

    # Save the analysis
    io.write_h5mu(mdata, "BackhedF_2015")


if __name__ == "__main__":
    run()
