from pathlib import Path
import pandas as pd


with open(Path(__file__).parent / "README.md", "r") as fh:
    long_description = fh.read()

inventory = pd.read_csv(Path(__file__).parent / "inventory.csv")
