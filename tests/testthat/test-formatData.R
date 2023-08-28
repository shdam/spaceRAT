test_that("formatData works", {

    data("pData_dmap", package = "spaceRATScaffolds")

    expect_equal(formatPheno(pData_dmap), pData_dmap)

    pData_dmap2 <- data.frame(
        "rownames" = rownames(pData_dmap),
        "cell_types" = pData_dmap$cell_types)
    expect_equal(formatPheno(pData_dmap2), pData_dmap)

    pData_dmap$na <- rep(NA, nrow(pData_dmap))
    expect_message(formatPheno(pData_dmap, "na"), ".*in phenotype table contain NA.*")
})

test_that("formatData errors checks", {

    data("pData_dmap", package = "spaceRATScaffolds")

    expect_error(formatPheno(pData_dmap, colname = "Wrong"), "Column Wrong was not found in annotation data.")
})

