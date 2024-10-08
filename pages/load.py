import streamlit as st
from mudata_explorer.helpers import save_load
from mudata_explorer.helpers import cirro
from mudata_explorer.app.sidebar import setup_sidebar


def run():
    setup_sidebar("load")

    st.subheader("Load Source Data + Figures")

    with st.container(border=1):
        st.write("#### Load From Cirro")
        cirro.load_from_cirro(
            filter_process_ids=["mudata-h5mu"]
        )

    save_load.upload_mdata(
        st.container(border=True)
    )

    # Let the user build a dataset by uploading tables
    with st.container(border=1):
        st.write("#### Load from Tables")
        st.page_link(
            "pages/load_tables.py",
            label="Build a MuData from Spreadsheets",
            icon=":material/upload:"
        )

    # Let the user browse existing datasets
    with st.container(border=1):
        st.write("#### Build Report")
        st.page_link(
            "pages/load_microbiome.py",
            label="Microbiome Report",
            icon=":material/microbiology:"
        )

    # Let the user load from a URL
    with st.container(border=1):
        st.write("#### Load from URL")
        url = st.text_input(
            "MuData URL (*.h5mu)",
            key="url",
            placeholder="https://example.com/data.h5mu"
        )
        if url:
            st.markdown(f"[Permalink](?file={url})")
            if st.button("Load"):
                st.session_state["file"] = url
                st.rerun()


if __name__ == "__main__":
    run()
