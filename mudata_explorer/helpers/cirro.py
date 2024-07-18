from typing import Optional, Union
from cirro import CirroApi, DataPortal, DataPortalProject
from cirro import DataPortalDataset
from cirro.auth.device_code import DeviceCodeAuth
from cirro.config import AppConfig, list_tenants
from io import StringIO
from muon import MuData
from mudata_explorer import app
from mudata_explorer.helpers.cirro_readers import util, mudata
import streamlit as st
from streamlit.delta_generator import DeltaGenerator
from streamlit.runtime.scriptrunner import script_run_context
from tempfile import TemporaryDirectory
from threading import Thread
from time import sleep


def setup_cirro_client():

    tenant_dict = {
        tenant['displayName']: tenant['domain']
        for tenant in list_tenants()
    }

    # Let the user select a tenant
    tenant = st.selectbox(
        "Organization",
        ["< select for login >"] + list(tenant_dict.keys())
    )

    domain = tenant_dict.get(tenant)

    if domain:
        _cirro_login(domain, st.empty())


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


def load_from_cirro():

    st.write("#### Load from Cirro")

    # If the Cirro client has not been set up,
    # prompt the user to set it up
    if not st.session_state.get("Cirro"):
        setup_cirro_client()
        return

    # Ask the user to select a project
    project = _select_project("load")
    if project is None:
        return

    # Ask the user to select a dataset
    dataset = _select_dataset(project)
    if dataset is None:
        return

    # Read the MuData object from the contents of the dataset
    mdata = _read_dataset(dataset)

    # If no data was read, stop
    if mdata is None:
        return

    # Get the hash of the data
    _, hash, _ = app.get_dat_hash(mdata)

    if st.button("Load Dataset"):
        app.hydrate_uns(mdata)
        app.set_mdata(mdata)
        app.set_mdata_hash(hash)
        st.page_link(
            "pages/views.py",
            label="View Data",
            icon=":material/insert_chart:"
        )


def save_to_cirro():

    st.write("#### Save to Cirro")

    # If the Cirro client has not been set up,
    # prompt the user to set it up
    if not st.session_state.get("Cirro"):
        setup_cirro_client()
        return

    # Ask the user to select a project
    project = _select_project("save")
    if project is None:
        return

    # Ask the user to provide a name and description for the dataset
    name = st.text_input("Name")
    description = st.text_area("Description")

    if not st.button("Save"):
        return

    if not name:
        st.error("Name is required")
        return

    if not description:
        st.error("Description is required")
        return

    # Get the active dataset
    mdata = app.get_mdata()

    # If there is no active dataset
    if mdata is None:
        st.error("No dataset to save")
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
            st.error(f"Error: {e}")
            return

    st.write(
        f"Wrote {name} to {project.name} ({size})"
    )


def _select_project(key: str) -> DataPortalProject:
    """Get the list of projects available from Cirro."""

    portal: DataPortal = st.session_state.get("Cirro")
    if not portal:
        util.clear_cirro_client()

    # Try to get the list of projects
    try:
        projects = portal.list_projects()
    # If there is an error
    except Exception as e:
        # Report it to the user
        st.error(f"Error: {e}")
        util.clear_cirro_client()

    # Give the user the option to select a project
    project = st.selectbox(
        "Project",
        ["< select a project >"] + [p.name for p in projects],
        key=f"{key}_project"
    )

    if project == "< select a project >":
        return

    # Return the project object
    return next(p for p in projects if p.name == project)


def _select_dataset(project: DataPortalProject) -> DataPortalDataset:

    # Try to get the list of datasets in the project
    try:
        datasets = project.list_datasets()
    except Exception as e:
        st.error(f"Error: {e}")
        util.clear_cirro_client()

    # Filter down to the datasets which have parsing configured
    datasets = [
        d for d in datasets
        if _read_dataset(d, check_only=True)
    ]

    # Give the user the option to select a dataset
    dataset = st.selectbox(
        "Dataset",
        ["< select a dataset >"] + [ds.name for ds in datasets],
        key="select-dataset"
    )

    if dataset == "< select a dataset >":
        return

    # Return the dataset object
    return next(ds for ds in datasets if ds.name == dataset)


def _read_dataset(
    dataset: DataPortalDataset,
    check_only: bool = False
) -> Union[bool, Optional[MuData]]:
    """
    Read a dataset from Cirro.
    Report any errors that arise.
    If check_only is True, only check if the dataset's type is valid
    """

    # MuData datasets
    if dataset.process_id == "mudata-h5mu":
        if check_only:
            return True
        else:
            return mudata.read(dataset)


# def autoretry(func, retries=15, exception=Exception):
#     def inner(*args, **kwargs):

#         result = None
#         for i in range(retries):
#             try:
#                 result = func(*args, **kwargs)
#             except exception as e:
#                 if i == (retries - 1):
#                     raise e
#                 else:
#                     sleep(0.1)
#             if result is not None:
#                 break
#         return result

#     return inner
