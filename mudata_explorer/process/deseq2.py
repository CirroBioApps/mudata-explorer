from typing import Optional
import numpy as np
import pandas as pd
from mudata_explorer.base.process import Process
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats


class RunDESeq2(Process):

    type = "deseq2"
    name = "DESeq2 Test"
    help_text = """
The DESeq2 test is a method for differential expression analysis of count data,
using shrinkage estimation for dispersions and fold changes to improve stability
and interpretability of estimates.
The Python implementation of DESeq2 is used here (`PyDESeq2`).

    """ # noqa

    category = "Statistical Tests"
    schema = {
        "table": {
            "type": "object",
            "label": "Data Table",
            "properties": {
                "data": {
                    "type": "dataframe",
                    "label": "Data",
                    "help": "Select the data which will be used to perform the Kruskal-Wallis test", # noqa
                    "select_columns": True,
                    "query": "",
                },
                "grouping": {
                    "type": "dataframe",
                    "label": "Grouping",
                    "help": "Select the column containing categories to compare", # noqa
                    "columns": {"grouping": {"label": "Grouping"}}
                }
            }
        },
        "comparison": {
            "type": "object",
            "label": "Comparison Between Groups",
            "properties": {
                "ref_level": {
                    "type": "string",
                    "label": "Reference Level",
                    "help": "Select the reference level for the comparison"
                },
                "comp_level": {
                    "type": "string",
                    "label": "Comparison Level",
                    "help": "Select the comparison level for the comparison"
                },
                "alpha": {
                    "type": "float",
                    "default": 0.05,
                    "label": "Alpha"
                },
                "independent_filter": {
                    "type": "boolean",
                    "default": True,
                    "label": "Independent Filtering"
                }
            }
        },
        "outputs": {
            "type": "object",
            "label": "Outputs",
            "properties": {
                "dest_key": {
                    "type": "string",
                    "default": "deseq2",
                    "label": "Label to use for results",
                    "help": """
                    Key to use when saving the output
                    """
                }
            }
        }
    }
    outputs = {
        "results": {
            "type": pd.DataFrame,
            "label": "DESeq2 Results",
            "desc": "Results of the DESeq2 Analysis.",
            "modality": "table.data.tables",
            "axis": "table.data.axis.T",
            "attr": "outputs.dest_key"
        }
    }

    def execute(self):

        df: Optional[pd.DataFrame] = self.params["table.data.dataframe"]
        if df is None:
            raise ValueError("Data table is empty.")
        df = df.dropna()
        assert df.shape[0] > 0, "Null values in all rows."

        group: pd.Series = (
            self.params
            ["table.grouping.dataframe"]
            ["grouping"]
            .dropna()
        )

        # Make sure that the comp_level and ref_level are in the grouping
        comp_level = self.params["comparison.comp_level"]
        assert comp_level and comp_level in group.unique(), \
            "Comparison level not found in grouping"
        ref_level = self.params["comparison.ref_level"]
        assert ref_level and ref_level in group.unique(), \
            "Reference level not found in grouping"

        # Get the set of observations which have both data and grouping
        shared = df.index.intersection(group.index)

        assert len(shared) > 0, "No shared rows between data and grouping"

        df = df.loc[shared]
        group = group.loc[shared]

        res = run_deseq2(
            df,
            group, 
            ref_level,
            comp_level,
            self.params["comparison.alpha"],
            self.params["comparison.independent_filter"]
        )

        self.save_results("results", res)


def run_deseq2(
    df: pd.DataFrame,
    group: pd.Series,
    ref_level: str,
    comp_level: str,
    alpha: float,
    independent_filter: bool
) -> pd.DataFrame:
    """
    Following the example shown:
    https://pydeseq2.readthedocs.io/en/latest/auto_examples/plot_step_by_step.html
    """

    assert comp_level and comp_level in group.unique(), \
        f"Comparison level not found in grouping ({comp_level})"
    assert ref_level and ref_level in group.unique(), \
        f"Reference level not found in grouping ({ref_level})"
    
    # Only keep those rows which are involved in the comparison
    group = group[group.isin([ref_level, comp_level])]
    df = df.loc[group.index]

    # Initialize the dataset
    dds = DeseqDataSet(
        counts=df.clip(lower=0).applymap(int),
        metadata=pd.DataFrame(dict(group=group)),
        design_factors="group",
        ref_level=["group", comp_level]
    )
    # Fit size factors
    dds.fit_size_factors()
    # Estimate dispersions
    dds.fit_genewise_dispersions()
    # Fit dispersion trend coefficients
    dds.fit_dispersion_trend()
    # Fit dispersion priors
    dds.fit_dispersion_prior()
    # MAP dispersions
    dds.fit_MAP_dispersions()
    # Fit log fold changes
    dds.fit_LFC()
    # Fit Cooks
    dds.calculate_cooks()
    if dds.refit_cooks:
        # Replace outlier counts
        dds.refit()

    # Run the DESeq2 test
    stat_res = DeseqStats(dds, alpha=alpha, independent_filter=independent_filter)

    # Wald tests
    stat_res.run_wald_test()

    # Independent filtering or p-value adjustment
    if stat_res.independent_filter:
        stat_res._independent_filtering()
    else:
        stat_res._p_value_adjustment()

    # Create a summary of the results
    stat_res.summary()

    # Calculate the -log10 p-value
    min_p = (
        stat_res.results_df
        .query("pvalue > 0")
        ["pvalue"].min()
    )
    res = stat_res.results_df.assign(
        neg_log10_pvalue=(
            -stat_res
            .results_df["pvalue"]
            .clip(lower=min_p)
            .apply(np.log10)
        )
    )

    return res