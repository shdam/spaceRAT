utils::data("DMAP_exprs", "DMAP_pData", package = "spaceRATScaffolds")
scaffold <- buildScaffold(DMAP_exprs, DMAP_pData, "cell_types", data = "logged")

test_that("loadingPlot works", {
    l <- loadingPlot(scaffold)
    df <- loadingPlot(scaffold, df_only = TRUE)
    expect_s3_class(l, "ggplot")
    expect_s3_class(df, "data.frame")
    expect_equal(nrow(df), 12)
    expect_equal(ncol(df), 4)
    expect_equal(colnames(df), c("PC1", "PC2", "class", "gene"))

    expect_s3_class(loadingPlot(space, dims = c(2,3)), "ggplot")
    expect_s3_class(loadingPlot(space, dims = c(2,3), num_genes = 5), "ggplot")
    expect_warning(loadingPlot(space, dims = c(2,3), gene_name = "entrez"))
    expect_s3_class(loadingPlot(space, dims = c(2,3), angle = 45), "ggplot")
})


test_that("loadingPlot returns righ colnames", {
    df <- loadingPlot(space, dims = c(2,3), df_only = TRUE, num_genes = 4)
    expect_equal(colnames(df), c("PC2", "PC3", "class", "gene"))
    expect_equal(nrow(df), 16)
    expect_s3_class(df, "data.frame")
})





