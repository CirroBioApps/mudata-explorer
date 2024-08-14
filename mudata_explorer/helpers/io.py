import json
import numpy as np
import muon as mu
from tempfile import NamedTemporaryFile
from typing import Union


def jsonify(dat):
    if isinstance(dat, (list, np.ndarray)):
        return [jsonify(val) for val in dat]
    elif isinstance(dat, dict):
        return {kw: jsonify(val) for kw, val in dat.items()}
    else:
        return dat


def validate_json(dat):
    """Validate that an object can be serialized to JSON"""
    try:
        return json.loads(
            json.dumps(
                jsonify(dat),
                sort_keys=True
            )
        )
    except Exception as e:
        print(dat)
        raise ValueError(f"Could not serialize object to JSON: {e}")


def json_safe(obj: Union[str, dict]):
    if isinstance(obj, str):
        return json.loads(obj)
    else:
        return obj


def hydrate_uns(mdata: mu.MuData):
    prefix = "mudata-explorer-"
    for suffix in ["views", "history", "provenance", "settings"]:
        kw = f"{prefix}{suffix}"
        if kw in mdata.uns:
            if isinstance(mdata.uns[kw], str):
                mdata.uns[kw] = json.loads(mdata.uns[kw])


def dehydrate_uns(mdata: mu.MuData):
    prefix = "mudata-explorer-"
    for suffix in ["views", "history", "provenance", "settings"]:
        kw = f"{prefix}{suffix}"
        if kw in mdata.uns:
            if not isinstance(mdata.uns[kw], str):
                mdata.uns[kw] = json.dumps(mdata.uns[kw], sort_keys=True)


def mdata_to_binary(mdata: mu.MuData) -> bytes:
    # Write out the MuData object to a temporary file
    with NamedTemporaryFile(suffix=".h5mu", delete=True) as tmp:

        # Write the MuData object to the file
        write_h5mu(mdata, tmp.file.name)

        # Get the file object in bytes
        return tmp.read()


def write_h5mu(mdata: mu.MuData, file_name: str):
    # Convert any .uns objects to strings
    dehydrate_uns(mdata)

    # Format as h5mu
    mdata.write(file_name)

    # Convert back to objects
    hydrate_uns(mdata)


def read_h5mu(h5mu_file) -> mu.MuData:

    with NamedTemporaryFile(suffix=".h5mu", delete=True) as tmp:
        with open(tmp.file.name, "wb") as f:
            f.write(h5mu_file.getvalue())
        mdata = mu.read_h5mu(tmp.file.name)

    hydrate_uns(mdata)

    return mdata


