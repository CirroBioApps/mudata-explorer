from mudata_explorer.public_data.curatedMetagenomicData.inventory import inventory
import streamlit as st


def run():
    st.markdown("""#### Gallery: curatedMetagenomicData

The source data for all of the Microbiome Reports in this gallery 
were obtained from the [curatedMetagenomicData project](https://waldronlab.io/curatedMetagenomicData/).
""")

    for _, r in inventory.iterrows():
        with st.container(border=True):
            cols = st.columns([2, 1])
            cols[0].write(f"""**{r['Dataset Name']}** ({r['Total Samples']})
                     
Compared By: **{r['Comparison By']}**
""")
            cols[1].write(f"![Thumbnail]({r['png']})")

            if cols[0].button("View", key=f"view-{r['Dataset Name']}"):
                st.session_state["file"] = r["path"]
                st.switch_page("pages/view_all.py")
