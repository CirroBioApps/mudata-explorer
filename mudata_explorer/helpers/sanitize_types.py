import pandas as pd
from streamlit.delta_generator import DeltaGenerator


def sanitize_types(df: pd.DataFrame, container: DeltaGenerator):
    # Convert every column to the type of the first non-null value
    to_drop = []
    for cname, col in df.items():
        if col.dtype == "object":
            try:
                df[cname] = col.astype(type(col.dropna().values[0]))
            except Exception:
                to_drop.append(cname)
        elif col.dtype == "str":
            to_drop.append(cname)
    if len(to_drop) > 0:
        container.write(f"Dropping non-numeric columns: {to_drop}")
        df = df.drop(columns=to_drop)

    return df

