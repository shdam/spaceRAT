utils::data("exprs_dmap", "pData_dmap", package = "spaceRATScaffolds")
space <- buildScaffold(exprs_dmap, pData_dmap, "cell_types", data = "logged")

test_that("loadingPlot works", {
    l <- loadingPlot(space)
    df <- loadingPlot(space, df_only = TRUE)
    expect_s3_class(l, "ggplot")
    expect_s3_class(df, "data.frame")
    expect_equal(nrow(df), 12)
    expect_equal(ncol(df), 4)
    expect_equal(colnames(df), c("PC1", "PC2", "class", "gene"))
})


test_that("loadingPlot returns righ colnames", {
    df <- loadingPlot(space, dims = c(2,3), df_only = TRUE, num_genes = 4)
    expect_equal(colnames(df), c("PC2", "PC3", "class", "gene"))
    expect_equal(nrow(df), 16)
    expect_s3_class(df, "data.frame")
})





