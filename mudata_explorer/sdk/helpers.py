def collapse_params(params: dict):
    """
    Params may be provided as nested dictionaries, but we want to
    flatten them to a single level. This function does that, using
    underscores to separate keys.
    """
    output = dict()
    for key, val in params.items():
        if isinstance(val, dict):
            for sub_key, sub_val in val.items():
                output[f"{key}_{sub_key}"] = sub_val
        else:
            output[key] = val

    output = {
        kw.replace(".", "_"): val
        for kw, val in output.items()
    }
    if any([isinstance(val, dict) for val in output.values()]):
        return collapse_params(output)
    else:
        return output
