data("exprs_dmap")
data("pData_dmap")

se <- SummarizedExperiment::SummarizedExperiment(assays = S4Vectors::SimpleList("exprs" = exprs_dmap), colData = pData_dmap)
test_that("preprocess works", {


    preprocessed <- preprocess(exprs_dmap, pheno = pData_dmap, colname = "cell_types")

    expect_equal(preprocessed[[1]], exprs_dmap)
    expect_equal(preprocessed[[2]], pData_dmap$cell_types)

    preprocessed <- suppressWarnings(preprocess(exprs_dmap, pheno = pData_dmap, colname = "cell_types", annotation = "entrez"))

    expect_no_error(as.numeric(rownames(preprocessed[[1]])))


    preprocessed <- preprocess(exprs_dmap)
    expect_null(preprocessed[[2]])

    preprocessed <- preprocess(se, assay = "exprs")
    expect_equal(preprocessed[[1]], exprs_dmap)
    expect_equal(preprocessed[[2]], pData_dmap$cell_types)


    preprocessed <- preprocess(se, assay = "exprs", classes = c("HSC", "ERY"))

    expect_equal(as.character(unique(preprocessed[[2]])), c("HSC", "ERY"))
    expect_equal(ncol(preprocessed[[1]]), 47)

})

test_that("warnings",{
    expect_error(preprocess(se, colname = "wrone_col"), "Column wrone_col wasn't found in pheno data.")
    expect_warning(preprocess(se, annotation = NULL), "Gene identifier has not been resolved.*")
})
