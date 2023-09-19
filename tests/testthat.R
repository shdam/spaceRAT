library(testthat)
library(spaceRAT)
library(spaceRATScaffolds)

utils::data(list = c("DMAP_exprs", "DMAP_pData", "gene_id_converter_hs", "ilaria_counts", "ilaria_pData"), package = "spaceRATScaffolds")

test_check("spaceRAT")
