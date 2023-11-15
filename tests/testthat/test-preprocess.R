data("DMAP_exprs", package = "spaceRATScaffolds")
data("DMAP_pData", package = "spaceRATScaffolds")

se <- SummarizedExperiment::SummarizedExperiment(assays = S4Vectors::SimpleList("exprs" = DMAP_exprs), colData = DMAP_pData)
test_that("preprocess works", {


    preprocessed <- preprocess(DMAP_exprs, data = "exprs", pheno = DMAP_pData, colname = "cell_types", threshold = 10)

    expect_equal(preprocessed[[1]], DMAP_exprs)
    expect_equal(preprocessed[[2]], DMAP_pData)

    preprocessed <- suppressWarnings(preprocess(DMAP_exprs, pheno = DMAP_pData, data = "exprs", colname = "cell_types", annotation = "entrez"))

    expect_no_error(as.numeric(rownames(preprocessed[[1]])))


    preprocessed <- preprocess(DMAP_exprs, data = "exprs")
    expect_null(preprocessed[[2]])

    preprocessed <- preprocess(se, assay = "exprs", data = "exprs")
    expect_equal(preprocessed[[1]], DMAP_exprs)
    expect_equal(preprocessed[[2]], DMAP_pData)


    preprocessed <- preprocess(se, assay = "exprs", classes = c("HSC", "ERY"), data = "exprs")

    expect_equal(as.character(unique(preprocessed[[2]]$cell_types)), c("HSC", "ERY"))
    expect_equal(ncol(preprocessed[[1]]), 47)

})

test_that("warnings",{
    expect_error(preprocess(se, assay = "exprs", colname = "wrone_col", data = "exprs"), "Column wrone_col was not found in pheno/colData data.")
    expect_warning(preprocess(se, assay = "exprs", annotation = NULL, data = "exprs"), "Gene identifier has not been resolved.*")
})
