library(testthat)
library(spaceRAT)

utils::data("exprs_dmap", "pData_dmap", "gene_id_converter_hs", "exprs_ilaria", "pData_ilaria", package = "spaceRAT")

test_check("spaceRAT")
