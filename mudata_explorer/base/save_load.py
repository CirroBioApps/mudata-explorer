import hashlib
import muon as mu
from mudata_explorer.base.base import MuDataAppHelpers
from streamlit.delta_generator import DeltaGenerator
from tempfile import NamedTemporaryFile
import json


class SaveLoad(MuDataAppHelpers):

    on_change: callable

    def __init__(self, on_change):
        self.on_change = on_change

    @property
    def key_prefix(self):
        return f"save-load-{self.refresh_ix('save-load')}-"

    def param_key(self, key: str):
        return f"{self.key_prefix}{key}"

    def show(self, empty: DeltaGenerator):
        container = empty.container()
        self.upload_button(container)
        self.download_button(container)

    def upload_button(self, container: DeltaGenerator):
        h5mu_file = container.file_uploader(
            "Upload MuData (.h5mu)",
            key=self.param_key("upload-h5mu")
        )
        if h5mu_file is None:
            return
        try:
            mdata = self.read_h5mu(h5mu_file)
        except ValueError as e:
            container.write(
                f"Could not parse file: {h5mu_file.name}\n\n{str(e)}"
            )
            return

        self.set_mdata(mdata)
        self.on_change()

    def read_h5mu(self, h5mu_file):
        with open(h5mu_file.name, "wb") as f:
            f.write(h5mu_file.getvalue())
        return mu.read(h5mu_file.name)

    def download_button(self, container: DeltaGenerator):

        self.summarize_mdata(container)

        mdata = self.get_mdata()
        if mdata is None:
            return

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

            container.write(f"Unique Hash: {hash}")

            # Name the downloaded file for the hash of the data
            container.download_button(
                "Download MuData",
                dat,
                f"mudata-{hash}.h5mu",
                help="Click here to download the MuData object as a file."
            )

    def hash_dat(self, dat):
        """Compute the hash of an object."""
        hash = hashlib.sha256()
        hash.update(dat)
        return hash.hexdigest()
