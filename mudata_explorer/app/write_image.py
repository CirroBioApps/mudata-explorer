import streamlit as st
import os
from tempfile import TemporaryDirectory
from mudata_explorer.base.view import View
from plotly import io


@st.dialog("Save Image", width='large')
def dialog_write_image(view: View):
    # Ask the user for the filename
    file_prefix = st.text_input(
        "File Name",
        value="image"
    )

    # Get other optional metadata used to format the image
    width = st.number_input(
        "Width (px)",
        value=io.kaleido.scope.default_width,
        help="The width of the exported image in layout pixels. If the scale property is 1.0, this will also be the width of the exported image in physical pixels."
    )

    height = st.number_input(
        "Height (px)",
        value=io.kaleido.scope.default_height,
        help="The height of the exported image in layout pixels. If the scale property is 1.0, this will also be the height of the exported image in physical pixels."
    )

    # Ask the user for the file type
    file_type = st.selectbox(
        "File Type",
        ["png", "jpeg", "svg", "pdf", "webp"],
        index=None
    )
    if file_type is None:
        return
    
    filename = f"{file_prefix}.{file_type}"
    with TemporaryDirectory() as tmp_dir:
        file_path = os.path.join(tmp_dir, filename)
        try:
            view.write_image(file_path, width=width, height=height)
        except Exception as e:
            st.exception(e)
        with open(file_path, "rb") as f:
            img = f.read()

    st.download_button(
        label=f"Download Image ({filename})",
        data=img,
        file_name=filename,
        mime=f"image/{file_type}"
    )

