data("exprs_dmap")
data("pData_dmap")

se <- SummarizedExperiment::SummarizedExperiment(assays = S4Vectors::SimpleList("exprs" = exprs_dmap), colData = pData_dmap)
test_that("preprocess works", {


    preprocessed <- preprocess(exprs_dmap, data = "logged", pheno = pData_dmap, colname = "cell_types", threshold = 10)

    expect_equal(preprocessed[[1]], exprs_dmap)
    expect_equal(preprocessed[[2]], pData_dmap)

    preprocessed <- suppressWarnings(preprocess(exprs_dmap, pheno = pData_dmap, data = "logged", colname = "cell_types", annotation = "entrez"))

    expect_no_error(as.numeric(rownames(preprocessed[[1]])))


    preprocessed <- preprocess(exprs_dmap, data = "logged")
    expect_null(preprocessed[[2]])

    preprocessed <- preprocess(se, assay = "exprs", data = "logged")
    expect_equal(preprocessed[[1]], exprs_dmap)
    expect_equal(preprocessed[[2]], pData_dmap)


    preprocessed <- preprocess(se, assay = "exprs", classes = c("HSC", "ERY"), data = "logged")

    expect_equal(as.character(unique(preprocessed[[2]]$cell_types)), c("HSC", "ERY"))
    expect_equal(ncol(preprocessed[[1]]), 47)

})

test_that("warnings",{
    expect_error(preprocess(se, assay = "exprs", colname = "wrone_col", data = "logged"), "Column wrone_col was not found in pheno/colData data.")
    expect_warning(preprocess(se, assay = "exprs", annotation = NULL, data = "logged"), "Gene identifier has not been resolved.*")
})
