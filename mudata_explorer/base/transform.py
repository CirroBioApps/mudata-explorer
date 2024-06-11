import numpy as np
import pandas as pd
from scipy import stats


class Transform:
    id: str
    name: str

    def run(df: pd.DataFrame) -> pd.DataFrame:
        pass


class FillNaZero(Transform):
    id = "fill_na_zero"
    name = "Fill Missing Values with 0"

    def run(df: pd.DataFrame) -> pd.DataFrame:
        return df.fillna(0)


class AddOne(Transform):
    id = "add_one"
    name = "Add 1 to All Values"

    def run(df: pd.DataFrame) -> pd.DataFrame:
        return df.apply(lambda x: x + 1)


class ZscoreRows(Transform):
    id = "zscores_rows"
    name = "Calculate Z-Scores for Rows"

    def run(df: pd.DataFrame) -> pd.DataFrame:
        return df.apply(safe_zscore, axis=1)


class ZscoreCols(Transform):
    id = "zscores_cols"
    name = "Calculate Z-Scores for Columns"

    def run(df: pd.DataFrame) -> pd.DataFrame:
        return df.apply(safe_zscore, axis=0)


def safe_zscore(x: pd.Series):
    """Calculate the z-score for a series, ignoring NaN values."""
    return (x - x.dropna().mean()) / x.dropna().std()


class LogTen(Transform):
    id = "log10"
    name = "Logarithm Base 10"

    def run(df: pd.DataFrame) -> pd.DataFrame:
        return df.apply(np.log10)


class LogTwo(Transform):
    id = "log2"
    name = "Logarithm Base 2"

    def run(df: pd.DataFrame) -> pd.DataFrame:
        return df.apply(np.log2)


class Log(Transform):
    id = "log"
    name = "Natural Logarithm"

    def run(df: pd.DataFrame) -> pd.DataFrame:
        return df.apply(np.log)


class Log1p(Transform):
    id = "log1p"
    name = "Natural Logarithm Plus 1"

    def run(df: pd.DataFrame) -> pd.DataFrame:
        return df.apply(np.log1p)


class Abs(Transform):
    id = "abs"
    name = "Absolute Value"

    def run(df: pd.DataFrame) -> pd.DataFrame:
        return df.apply(np.abs)


class Negative(Transform):
    id = "negative"
    name = "Negative Values (Multiply by -1)"

    def run(df: pd.DataFrame) -> pd.DataFrame:
        return df.apply(lambda v: v * -1)
