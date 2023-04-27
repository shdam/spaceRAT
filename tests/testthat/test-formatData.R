test_that("formatData works", {

    # data("pData_dmap")

    expect_equal(formatPheno(pData_dmap), pData_dmap)

    expect_error(formatPheno(pData_dmap, colname = "Wrong"), "Column Wrong was not found in annotation data.")

    pData_dmap$na <- rep(NA, nrow(pData_dmap))
    expect_message(formatPheno(pData_dmap, "na"), ".*in phenotype table contain NA.*")
})
