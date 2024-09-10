import pandas as pd
from pathlib import Path
import json
from os.path import getctime
import click

repo = "https://github.com/CirroBioApps/mudata-examples/raw/main"


def describe(config: dict, ix: int, basename: str, examples_repo: str) -> dict:

    df = pd.read_csv(
        basename + ".relative_abundance.tsv",
        sep="\t",
        index_col=0
    )

    return {
        "Dataset Name": config["dataset_name"],
        "Total Samples": n_samples(config, df),
        "Comparison By": metadata(config, df),
        "Dataset": find_file(basename, ix, examples_repo),
        "n": df.shape[0]
    }


def find_file(basename: str, ix: int, examples_repo: str) -> str:
    folder, prefix = basename.rsplit("/", 1)
    files = list(Path(folder).rglob(f"{prefix}-{ix}*.h5mu"))
    if len(files) == 0:
        return
    assert len(files) > 0, f"No files found for {basename}-{ix}"
    # Use the newest file
    latest_file = max(files, key=getctime)
    rel_path = Path(latest_file).relative_to(Path(examples_repo))
    path = f"{repo}/{rel_path}"

    return f"[**{latest_file.name}**](https://mudata-explorer.streamlit.app/views?file={path})"


def n_samples(config: dict, df: pd.DataFrame) -> str:
    ntot = df.shape[0]
    if 'query' in config:
        nfilt = df.query(config['query']).shape[0]
        return f"{nfilt:,} of {ntot:,} samples where {config['query']}"
    else:
        return f"{ntot:,} samples"


def metadata(config: dict, df: pd.DataFrame) -> str:
    assert config["compare_by"] in df.columns, \
        f"Column {config['compare_by']} not found in {df.columns}"
    vals = df[config["compare_by"]].dropna()

    if config.get('is_categorical'):
        vc = vals.value_counts()

        if vc.shape[0] > 5:
            n_other = vc.iloc[5:].sum()
            vc = vc.iloc[:5]
            vc["Other"] = n_other

        counts = ", ".join(
            [f"{k}: {v:,}" for k, v in vc.items()]
        )
        return f"{config['label']} - {counts}"
    else:
        return config['label']


def write(inventory: pd.DataFrame):
    with open("README.md", "w") as f:

        f.write("""
The curatedMetagenomicData package provides standardized, curated human microbiome data for novel analyses.
For complete information on this project, 
[see their official source material](https://waldronlab.io/curatedMetagenomicData/articles/curatedMetagenomicData.html).

These datasets have been used to demonstrate the functionality of the MuData Explorer.
A standard set of analyses have been performed on each dataset, each comparing the microbiome
samples based on a metadata category provided by the authors.

![HMP IBD 2019 UMAP](https://github.com/CirroBioApps/mudata-explorer/raw/main/mudata_explorer/public_data/curatedMetagenomicData/screenshots/HMP_2019_ibdmdb-0-study_condition-faee1afa1755d7ba.UMAP.png)

DISCLAIMER: The datasets are provided as-is, and the analyses are for demonstration purposes only.
None of the results should be considered scientifically valid or accurate, and are not a representation
of the original authors' work, or the conclusions of the curatedMetagenomicData project.

""")
        f.write("## Datasets\n\n")
        inventory.to_markdown(f, index=False)


@click.command()
@click.argument('examples_repo')
def run(examples_repo):
    inventory = (
        pd.DataFrame([
            describe(config, ix, str(config_file).replace(".config.json", ""), examples_repo)
            for config_file in (Path(examples_repo) / "curatedMetagenomicData/data").rglob("*.config.json")
            for ix, config in enumerate(json.load(config_file.open()))
        ])
        .sort_values(by="n", ascending=False)
        .drop(columns=["n"])
    )

    # Write out as markdown
    write(inventory)


if __name__ == "__main__":
    run()
