test_that("findDEGenes() returns a vector of gene names",{
    object <- spaceRAT:::preprocess(exprs_dmap,
                                  colname = "cell_types",
                                  pheno = pData_dmap)

    counts_scaffold <- object[[1]]
    cell_types <- object[[2]]
    rm(object)
    DEgenes <- spaceRAT:::findDEGenes(counts_scaffold, cell_types, 0.05,2)
    expect_type(DEgenes,"character")
})
