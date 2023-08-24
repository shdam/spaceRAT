library(testthat)
library(spaceRAT)
library(spaceRATScaffolds)

utils::data(list = c("exprs_dmap", "pData_dmap", "gene_id_converter_hs", "counts_ilaria", "pData_ilaria"), package = "spaceRATScaffolds")

test_check("spaceRAT")
