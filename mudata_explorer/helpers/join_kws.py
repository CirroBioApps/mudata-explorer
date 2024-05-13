def join_kws(*kws):
    for kw in kws:
        if kw is not None:
            msg = f"Expected string, got {type(kw)} ({kw})."
            assert isinstance(kw, str), msg
    return '.'.join([kw for kw in kws if kw is not None])
