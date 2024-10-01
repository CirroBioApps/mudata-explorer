import json
import muon as mu
from pathlib import Path

# Load the explanatory dataset
explanation = mu.read_h5mu(Path(__file__).parent / "microbiome-report-176ff7077db10524.h5mu")

for suffix in ["history", "views", "provenence", "settings", "history"]:
    kw = f"mudata-explorer-{suffix}"
    if kw in explanation.uns:
        if isinstance(explanation.uns[kw], str):
            explanation.uns[kw] = json.loads(explanation.uns[kw])
