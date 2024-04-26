import hashlib
import mudata as mu
from mudata_explorer.base.view import View
from streamlit.delta_generator import DeltaGenerator
from tempfile import NamedTemporaryFile
import json


class DownloadMuData(View):

    type = "download-mudata"
    name = "Download MuData"
    desc = "Save the MuData as a local file."
    categories = ["Data Processing"]
    defaults = {}

    def display(self, container: DeltaGenerator):

        mdata = self.get_mdata()

        # Convert the view data to a JSON string
        mdata.uns["mudata-explorer-views"] = json.dumps(
            mdata.uns["mudata-explorer-views"]
        )

        # Write out the MuData object to a temporary file
        with NamedTemporaryFile(suffix=".h5mu", delete=True) as tmp:
            # Format as h5mu
            mdata.write(tmp.file.name)

            # Get the file object in bytes
            dat = tmp.read()

            # Compute the hash of the file
            hash = self.hash_dat(dat)

            # Name the downloaded file for the hash of the data
            container.download_button(
                "Download MuData",
                dat,
                f"mudata-{hash}.h5mu",
                help="Click here to download the MuData object as a file."
            )
            container.write(tmp.name)

    def hash_dat(self, dat):
        """Compute the hash of an object."""
        hash = hashlib.sha256()
        hash.update(dat)
        return hash.hexdigest()
