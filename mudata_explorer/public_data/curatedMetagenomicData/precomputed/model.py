from mudata_explorer.public_data.curatedMetagenomicData.inventory import inventory
from mudata_explorer.app.url import load_url
import streamlit as st


def run():

    st.markdown("""**curatedMetagenomicData**

The _curatedMetagenomicData_ package provides standardized, curated human microbiome data for novel analyses.
It includes gene families, marker abundance, marker presence, pathway abundance, pathway coverage, and
relative abundance for samples collected from different body sites. The bacterial, fungal, and archaeal
taxonomic abundances for each sample were calculated with MetaPhlAn3, and metabolic functional potential
was calculated with HUMAnN3.

For complete information on this project,
[click here](https://waldronlab.io/curatedMetagenomicData/).
                
A list of all available studies can be found here: [Available Studies](https://waldronlab.io/curatedMetagenomicData/articles/available-studies.html)
                """)

    # Format a nice display name for the inventory
    inventory_dict = {
        f"{r['Dataset Name']} ({r['Total Samples']}): Comparison By: {r['Comparison By']}": r["path"]
        for _, r in inventory.iterrows()
    }

    # Let the user pick one of the precomputed inventories
    selection = st.selectbox(
        "Load a pre-analyzed microbiome dataset",
        options=inventory_dict,
        index=None,
        placeholder="Select a dataset to load"
    )

    if selection is not None:
        url = inventory_dict[selection]
        load_url(url)
        st.switch_page("pages/view_all.py")

    