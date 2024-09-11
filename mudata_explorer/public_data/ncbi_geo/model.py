import anndata as ad
import GEOparse
import streamlit as st
from tempfile import TemporaryDirectory
from mudata_explorer.parsers.gene_expression import parse_adata
from mudata_explorer.apps.helpers import load_mdata


@st.cache_resource
def read_geo(geo_accession) -> ad.AnnData:

    # Load the GEO dataset
    with TemporaryDirectory() as tmp:
        gse = GEOparse.get_GEO(geo=geo_accession, destdir=tmp)

    # Get the expression data and set up the unique gene IDs as the index
    expression_df = (
        gse.table
        .assign(GENE_ID=gse.table.apply(
            lambda r: f"{r['IDENTIFIER']} ({r['ID_REF']})",
            axis=1
        ))
        .set_index('GENE_ID')
        .drop(columns=["IDENTIFIER", "ID_REF"])
    )

    # Create an AnnData object
    adata = ad.AnnData(
        X=expression_df.T,
        obs=gse.columns
    )

    return adata, gse.metadata


def run():
    st.markdown("""**NCBI GEO Datasets**

Read a dataset directly from the NCBI GEO database.

The full list of available datasets can be found in the [GEO Dataset Browser](https://www.ncbi.nlm.nih.gov/sites/GDSbrowser/).
""")

    geo_accession = st.text_input("Enter a GEO accession", placeholder="GDS****")

    if geo_accession is not None and len(geo_accession) > 0:

        if not geo_accession.startswith("GDS"):
            st.error("Invalid GEO accession. Please enter a valid GEO accession starting with 'GDS'.")
            return

        # Read the GEO dataset
        adata, metadata = read_geo(geo_accession)

        # Drop any genes which are all NaN
        adata = adata[:, adata.to_df().notnull().any(axis=0)]

        # Fill in any NaNs with 0
        adata.X = adata.to_df().fillna(0).values

        # Report the number of samples and features
        st.write(f"""**{metadata['title'][0]}**
                 
{metadata['description'][0]}

Read in data for **{adata.shape[0]:,} samples** and **{adata.shape[1]:,} genes**""")

        mdata = parse_adata(adata)
        if mdata is not None:
            load_mdata(mdata)
