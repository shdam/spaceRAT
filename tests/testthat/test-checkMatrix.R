test_that("checkMatrix returns an error when too many non-numerics are present", {
    mat <- data.frame(A = c("gene1", "gene2", "gene3"),
                      B = c(1, 2, "notnumeric"))
    expect_error(checkMatrix(mat), "Too many non-numerics in counts data! Please separate from counts.")
})

test_that("checkMatrix removes additional column when it is intended as rownames", {
    mat <- data.frame(gene = c("gene1", "gene2", "gene3"),
                      cond1 = c(1, 2, 3),
                      cond2 = c(3, 2, 1),
                      stringsAsFactors = FALSE)
    # colnames(mat)[1] <- "gene names"
    result <- checkMatrix(mat)
    expect_equal(rownames(result), c("gene1", "gene2", "gene3"))
    expect_equal(ncol(result), 2)
})


test_that("checkMatrix converts tibble to matrix", {
    mat <- tibble::tibble(gene = c("gene1", "gene2", "gene3"),
                          cond1 = c(1, 2, 3),
                          cond2 = c(3, 2, 1))
    result <- checkMatrix(mat)
    expect_true(is.matrix(result))
    expect_equal(row.names(result), c("gene1", "gene2", "gene3"))
})

test_that("checkMatrix handles matrix with rownames", {
    mat <- matrix(c(1, 2, 3, 3, 2, 1), ncol = 2)
    rownames(mat) <- c("gene1", "gene2", "gene3")
    result <- checkMatrix(mat)
    expect_identical(result, mat)
    expect_true(is.matrix(result))
    expect_equal(row.names(result), c("gene1", "gene2", "gene3"))
})

test_that("checkMatrix adds rownames", {
    mat <- matrix(c(1, 2, 3, 3, 2, 1), ncol = 2)
    expect_error(checkMatrix(mat), "The expression matrix does not contain rownames")
})
