import pandas as pd
from pathlib import Path
import json
from os.path import getctime

repo = "https://github.com/CirroBioApps/mudata-explorer/raw/main"


def describe(config: dict, ix: int, basename: str) -> dict:

    df = pd.read_csv(
        basename + ".relative_abundance.tsv",
        sep="\t",
        index_col=0
    )

    return {
        "Dataset Name": config["dataset_name"],
        "Total Samples": n_samples(config, df),
        "Comparison By": metadata(config, df),
        "Dataset": find_file(basename, ix),
        "n": df.shape[0]
    }


def find_file(basename: str, ix: int) -> str:
    folder, prefix = basename.rsplit("/", 1)
    files = list(Path(folder).rglob(f"{prefix}-{ix}*.h5mu"))
    if len(files) == 0:
        return
    assert len(files) > 0, f"No files found for {basename}-{ix}"
    # Use the newest file
    latest_file = max(files, key=getctime)
    rel_path = Path(latest_file).relative_to(Path(".."))
    path = f"{repo}/demo_data/curatedMetagenomicData/{rel_path}"

    return f"[**{latest_file.name}**](https://mudata-explorer.streamlit.app/views?file={path})"


def n_samples(config: dict, df: pd.DataFrame) -> str:
    ntot = df.shape[0]
    if 'query' in config:
        nfilt = df.query(config['query']).shape[0]
        return f"{nfilt:,} of {ntot:,} samples where {config['query']}"
    else:
        return f"{ntot:,} samples"


def metadata(config: dict, df: pd.DataFrame) -> str:
    assert config["cname"] in df.columns, \
        f"Column {config['cname']} not found in {df.columns}"
    vals = df[config["cname"]].dropna()

    if config.get('transform') == "quartiles":
        vals = pd.qcut(vals, q=4, labels=[f"Q{i+1}" for i in range(4)])

    vc = vals.value_counts()

    if vc.shape[0] > 5:
        n_other = vc.iloc[5:].sum()
        vc = vc.iloc[:5]
        vc["Other"] = n_other

    counts = ", ".join(
        [f"{k}: {v:,}" for k, v in vc.items()]
    )
    return f"{config['label']} - {counts}"


def write(inventory: pd.DataFrame):
    with open("README.md", "w") as f:
        f.write("# Example Datasets: curatedMetagenomicData\n\n")
        f.write("""
The curatedMetagenomicData package provides standardized, curated human microbiome data for novel analyses.
For complete information on this project, 
[see their official source material](https://waldronlab.io/curatedMetagenomicData/articles/curatedMetagenomicData.html).

These datasets have been used to demonstrate the functionality of the MuData Explorer.
A standard set of analyses have been performed on each dataset, each comparing the microbiome
samples based on a metadata category provided by the authors.

![HMP 2012 body site UMAP](https://github.com/CirroBioApps/mudata-explorer/raw/main/demo_data/curatedMetagenomicData/screenshots/HMP_2012-0-body_subsite-2f99b2563ed9e516.UMAP.png)

DISCLAIMER: The datasets are provided as-is, and the analyses are for demonstration purposes only.
None of the results should be considered scientifically valid or accurate, and are not a representation
of the original authors' work, or the conclusions of the curatedMetagenomicData project.

""")
        f.write("## Datasets\n\n")
        inventory.to_markdown(f, index=False)


if __name__ == "__main__":

    inventory = (
        pd.DataFrame([
            describe(config, ix, str(config_file).replace(".config.json", ""))
            for config_file in Path("../data").rglob("*.config.json")
            for ix, config in enumerate(json.load(config_file.open()))
        ])
        .sort_values(by="n", ascending=False)
        .drop(columns=["n"])
    )

    # Write out as markdown
    write(inventory)
