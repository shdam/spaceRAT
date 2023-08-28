# Create example data


# Test for logged data with some genes having total counts less than threshold
test_that("preFilter() removes genes with total counts less than threshold for logged data", {
    pre_filtered_counts <- preFilter(exprs_dmap, data = "logged", threshold = 100)
    expect_equal(dim(pre_filtered_counts), c(7940, 109))
    expect_false(any(rowSums(exp(pre_filtered_counts)) < 100))
})

# Test for logged data with all genes having total counts less than threshold
test_that("preFilter() throws an error when all genes have total counts less than threshold for logged data", {
    expect_error(preFilter(exprs_dmap, data = "logged", threshold = 100000),
                 "Low quality data! All genes have total counts less than *")
})

# Test for raw data with negative values
test_that("preFilter() throws an error for raw data with negative values", {
    expect_error(preFilter(exprs_dmap * -1, data = "counts", threshold = 10),
                 "Negative values are not allowed in raw count matrix!")
})

# Test for raw data with some genes having total counts less than threshold
test_that("preFilter() removes genes with total counts less than threshold for raw data", {
    pre_filtered_counts <- preFilter(exp(exprs_dmap), data = "counts", threshold = 100)
    expect_equal(dim(pre_filtered_counts), c(7940, 109))
    expect_false(any(rowSums(exp(pre_filtered_counts)) < 100))
})

# Test for raw data with all genes having total counts less than threshold
test_that("preFilter() throws an error when all genes have total counts less than threshold for raw data", {
    expect_error(preFilter(matrix(rep(0, 16), nrow = 4), data = "counts", threshold = 20),
                 "Low quality data! All genes have total counts less than*")
})

# Test for correct data input
test_that("preFilter() throws an error when all genes have total counts less than threshold for raw data", {
    expect_error(preFilter(matrix(rep(0, 16), nrow = 4), data = "whatnow", threshold = 20),
                 "Invalid 'data' argument. Please choose 'logged' or 'counts'.")
})
