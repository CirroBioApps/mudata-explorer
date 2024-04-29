import pandas as pd
from streamlit.delta_generator import DeltaGenerator


def sanitize_types(df: pd.DataFrame, container: DeltaGenerator, keep_str=False):
    # Convert every column to the type of the first non-null value
    to_drop = []
    for cname in df.columns:
        if df[cname].dtype == object:
            if isinstance(df[cname].dropna().values[0], str):
                if keep_str:
                    df[cname] = df[cname].fillna("").apply(str)
                    df.dtypes[cname] = str
                else:
                    to_drop.append(cname)
            else:
                try:
                    df[cname] = df[cname].apply(type(df[cname].dropna().values[0]))
                except Exception:
                    to_drop.append(cname)
        elif df[cname].dtype == str:
            to_drop.append(cname)
    if len(to_drop) > 0:
        container.write(f"Dropping non-numeric columns: {to_drop}")
        df = df.drop(columns=to_drop)

    return df
