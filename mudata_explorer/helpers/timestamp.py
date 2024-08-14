import pandas as pd


def get_timestamp():
    return str(pd.Timestamp.now())
