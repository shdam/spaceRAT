test_that("findDEGenes() returns a vector of gene names",{
  eset_dmap <- createEset(exprs_dmap,pData_dmap,"cell_types")
  DEgenes <- findDEGenes(eset_dmap,0.05,2)
  expect_type(DEgenes,"character")
})
