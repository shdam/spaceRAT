data("DMAP_pData", package = "spaceRATScaffolds")
test_that("formatData works", {
    expect_equal(formatPheno(DMAP_pData), DMAP_pData)

    DMAP_pData2 <- data.frame(
        "rownames" = rownames(DMAP_pData),
        "cell_types" = DMAP_pData$cell_types)
    expect_equal(formatPheno(DMAP_pData2), DMAP_pData)

    DMAP_pData$na <- rep(NA, nrow(DMAP_pData))
    expect_message(formatPheno(DMAP_pData, "na"), ".*in phenotype table contain NA.*")
})

test_that("formatData errors checks", {


    expect_error(formatPheno(DMAP_pData, colname = "Wrong"), "Column Wrong was not found in annotation data.")
})

