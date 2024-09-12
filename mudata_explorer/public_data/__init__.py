import mudata_explorer.public_data.curatedMetagenomicData.reanalyze.model as curatedMetagenomicData_reanalyze
import mudata_explorer.public_data.curatedMetagenomicData.precomputed.model as curatedMetagenomicData_precomputed
import mudata_explorer.public_data.ncbi_geo.model as ncbi_geo

repositories = {
    'NCBI GEO': ncbi_geo,
    'curatedMetagenomicData - Pre-analyzed': curatedMetagenomicData_precomputed,
    'curatedMetagenomicData - Run analysis': curatedMetagenomicData_reanalyze
}
