from mudata_explorer import app
import streamlit as st
from streamlit.delta_generator import DeltaGenerator


@st.experimental_dialog("Upload MuData")
def upload_button():
    h5mu_file = st.file_uploader(
        "**Upload MuData (.h5mu)**"
    )
    if h5mu_file is None:
        return
    try:
        mdata = app.read_h5mu(h5mu_file)
    except ValueError as e:
        st.write(
            f"Could not parse file: {h5mu_file.name}\n\n{str(e)}"
        )
        return

    if app.hash_dat(h5mu_file.read()) in h5mu_file.name:
        st.write("**Data Validated**: Unique hash matches file name.")
    else:
        st.write("Unique file hash not found in file name.")

    if st.button("Load Dataset"):
        app.set_mdata(mdata)
        st.rerun()


def show_hash(container: DeltaGenerator):

    # Get the current dataset along with its unique hash
    dat, hash, _ = app.get_dat_hash()
    if dat is None:
        return

    container.write(f"**Unique Hash**: {hash}")


@st.experimental_dialog("Download MuData")
def download_button():

    # Get the current dataset along with its unique hash
    dat, hash, size = app.get_dat_hash()
    if dat is None:
        return

    st.write(f"File is {size}")
    st.write(f"Unique hash: {hash}")

    # Name the downloaded file for the hash of the data
    if st.download_button(
        "Download MuData (.h5mu)",
        dat,
        f"mudata-{hash}.h5mu",
        help="Click here to download the MuData object as a file."
    ):
        st.rerun()
