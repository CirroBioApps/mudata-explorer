from pathlib import Path

with open(Path(__file__).parent / "README.md", "r") as fh:
    long_description = fh.read()
