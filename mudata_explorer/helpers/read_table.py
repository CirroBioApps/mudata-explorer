from streamlit.delta_generator import DeltaGenerator
import pandas as pd


def read_table(file, container: DeltaGenerator):

    if file is None:
        return

    if file.name.endswith("xlsx"):
        try:
            df = pd.read_excel(file)
        except pd.errors.EmptyDataError:
            container.write("Could not read data from file")
            return
    elif file.name.endswith("csv"):
        try:
            df = pd.read_csv(file)
        except pd.errors.EmptyDataError:
            container.write("Could not read data from file")
            return
    elif file.name.endswith("tsv"):
        try:
            df = pd.read_csv(file, sep="\t")
        except pd.errors.EmptyDataError:
            container.write("Could not read data from file")
            return
    else:
        container.error("File must be a CSV or TSV.")
        return

    # The first column must only have unique values
    if not df.iloc[:, 0].is_unique:
        container.error("The first column must have unique values.")
        return

    # Return the data
    return df
