
def nest_params(params: dict):
    output = dict()
    for key, value in params.items():
        if '.' in key:
            keys = key.split('.', 1)
            if keys[0] not in output:
                output[keys[0]] = dict()
            try:
                output[keys[0]][keys[1]] = value
            except TypeError as e:
                pass
            except Exception as e:
                print("Error nesting:", keys, value)
                raise e
        else:
            output[key] = value
    return {
        key: nest_params(value) if isinstance(value, dict) else value
        for key, value in output.items()
    }
