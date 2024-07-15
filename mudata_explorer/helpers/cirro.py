from cirro import CirroApi, DataPortal, DataPortalProject
from cirro import DataPortalDataset
from cirro.sdk.file import DataPortalFile
from cirro.auth.device_code import DeviceCodeAuth
from cirro.config import AppConfig, list_tenants
from io import StringIO
from muon import read_h5mu
from mudata_explorer import app
import streamlit as st
from streamlit.delta_generator import DeltaGenerator
from streamlit.runtime.scriptrunner import script_run_context
from tempfile import TemporaryDirectory
from threading import Thread
from time import sleep


def setup_cirro_client(container: DeltaGenerator):

    container.write("#### Cirro Login")

    tenant_dict = {
        tenant['displayName']: tenant['domain']
        for tenant in list_tenants()
    }

    # Let the user select a tenant
    tenant = container.selectbox(
        "Organization",
        ["< select for login >"] + list(tenant_dict.keys())
    )

    domain = tenant_dict.get(tenant)

    if domain:
        _cirro_login(domain, container.empty())


def _cirro_login(domain: str, container: DeltaGenerator):

    # Connect to Cirro - capturing the login URL
    auth_io = StringIO()
    cirro_login_thread = Thread(
        target=_cirro_login_sub,
        args=(auth_io, domain)
    )
    script_run_context.add_script_run_ctx(cirro_login_thread)

    cirro_login_thread.start()

    login_string = auth_io.getvalue()

    while len(login_string) == 0 and cirro_login_thread.is_alive():
        sleep(1)
        login_string = auth_io.getvalue()

    container.write(login_string)
    cirro_login_thread.join()
    container.empty()
    st.rerun()


def _cirro_login_sub(auth_io: StringIO, base_url: str):

    app_config = AppConfig(base_url=base_url)

    st.session_state['Cirro-auth_info'] = DeviceCodeAuth(
        region=app_config.region,
        client_id=app_config.client_id,
        auth_endpoint=app_config.auth_endpoint,
        enable_cache=False,
        auth_io=auth_io
    )

    st.session_state['Cirro-client'] = CirroApi(
        auth_info=st.session_state['Cirro-auth_info'],
        base_url=base_url
    )
    st.session_state['Cirro'] = DataPortal(
        client=st.session_state['Cirro-client']
    )


def _clear_cirro_client():
    # Wait for a second
    sleep(1)
    # Clear the client
    del st.session_state["Cirro"]
    # Let the user log in again
    st.rerun()


def load_from_cirro(container: DeltaGenerator):

    container.write("#### Load from Cirro")

    # Ask the user to select a project
    project = _select_project(container, "load")
    if project is None:
        return

    # Ask the user to select a dataset
    dataset = _select_dataset(container, project)
    if dataset is None:
        return

    # Pick the file from the dataset
    h5mu_file: DataPortalFile = _select_file(container, dataset)
    if h5mu_file is None:
        return

    # Download to a temporary folder
    with TemporaryDirectory() as tmp:
        # Download the file
        h5mu_file.download(tmp)
        # Path of the downloaded file
        filename = f"{tmp}/{h5mu_file.name}"
        # Read the MuData object
        mdata = read_h5mu(filename)
        # Calculate the hash of the object
        mdata_hash = app.hash_dat(open(filename, "rb").read())

    if mdata_hash in h5mu_file.name:
        container.write("**Data Validated**: Unique hash matches file name.")
    else:
        container.write("Unique file hash not found in file name.")

    if container.button("Load Dataset"):
        app.hydrate_uns(mdata)
        app.set_mdata(mdata)
        app.set_mdata_hash(mdata_hash)
        container.page_link("pages/views.py", label="View Data")


def save_to_cirro(container: DeltaGenerator):

    container.write("#### Save to Cirro")

    # Ask the user to select a project
    project = _select_project(container, "save")
    if project is None:
        return

    # Ask the user to provide a name and description for the dataset
    name = container.text_input("Name")
    description = container.text_area("Description")

    if not container.button("Save"):
        return

    if not name:
        container.error("Name is required")
        return

    if not description:
        container.error("Description is required")
        return

    # Get the active dataset
    mdata = app.get_mdata()

    # If there is no active dataset
    if mdata is None:
        container.error("No dataset to save")
        return

    # Get the binary blob, hash, and file size
    blob, hash, size = app.get_dat_hash(mdata)

    # Set the file name
    fn = f"{name.replace(' ', '-').lower()}-{hash}.h5mu"

    # Write out the dataset to a temporary file
    # and upload it to Cirro
    with TemporaryDirectory() as tmp:

        # Write the MuData object to the file
        open(f"{tmp}/{fn}", "wb").write(blob)

        # Upload the file to Cirro
        try:
            project.upload_dataset(
                name=name,
                description=description,
                process="mudata-h5mu",
                upload_folder=tmp
            )
        except Exception as e:
            container.error(f"Error: {e}")
            return

    container.write(
        f"Wrote {name} to {project.name} ({size})"
    )


def _select_project(container: DeltaGenerator, key: str) -> DataPortalProject:
    """Get the list of projects available from Cirro."""

    portal: DataPortal = st.session_state.get("Cirro")
    if not portal:
        _clear_cirro_client()

    # Try to get the list of projects
    try:
        projects = portal.list_projects()
    # If there is an error
    except Exception as e:
        # Report it to the user
        container.error(f"Error: {e}")
        _clear_cirro_client()

    # Give the user the option to select a project
    project = container.selectbox(
        "Project",
        ["< select a project >"] + [p.name for p in projects],
        key=f"{key}_project"
    )

    if project == "< select a project >":
        return

    # Return the project object
    return next(p for p in projects if p.name == project)


def _select_dataset(
    container: DeltaGenerator,
    project: DataPortalProject
) -> DataPortalDataset:

    # Try to get the list of datasets in the project
    try:
        datasets = project.list_datasets()
    except Exception as e:
        container.error(f"Error: {e}")
        _clear_cirro_client()

    # Filter down to the datasets which have parsing configured
    datasets = [
        d for d in datasets
        if _read_dataset(container, d, check_only=True)
    ]

    # Give the user the option to select a dataset
    dataset = container.selectbox(
        "Dataset",
        ["< select a dataset >"] + [ds.name for ds in datasets],
        key="select-dataset"
    )

    if dataset == "< select a dataset >":
        return

    # Return the dataset object
    return next(ds for ds in datasets if ds.name == dataset)


def _select_file(
    container: DeltaGenerator,
    dataset: DataPortalDataset
) -> DataPortalFile:

    # Try to get the list of files in the dataset
    try:
        files = dataset.list_files()
    except Exception as e:
        container.error(f"Error: {e}")
        _clear_cirro_client()

    # Give the user the option to select a file
    # (strip off the 'data/' prefix)
    file = container.selectbox(
        "File",
        ["< select a file >"] + [
            ds.name[5:] for ds in files
        ],
        key="select-file"
    )

    if file == "< select a file >":
        return

    # Return the file object
    return next(f for f in files if f.name[5:] == file)


def _read_dataset(
    container: DeltaGenerator,
    dataset: DataPortalDataset,
    check_only: bool = False
):
    """
    Read a dataset from Cirro.
    Report any errors that arise.
    If check_only is True, only check if the dataset's type is valid
    """

    # MuData datasets
    if dataset.process_id == "mudata-h5mu":
        if check_only:
            return True

    return False

