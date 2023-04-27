# Test case 1: Check output when some samples have no expression data
test_that("function removes samples with no expression data", {
    counts1 <- data.frame(A = c(1, 2, 3, 4, 5),
                          B = c(4, 5, 6, 7, 8), C = c(9, 10, 11, 12, 13))
    colnames(counts1) <- c("sample1", "sample5", "sample6")
    pheno1 <- data.frame(id = c("sample1", "sample2", "sample3"),
                         group = c("control", "treatment", "treatment"))
    rownames(pheno1) <- pheno1$id
    expected_output <- pheno1[-c(2:3), ,drop = FALSE]
    expect_message(
        missingExpression(pheno1, counts1))
    expect_equal(missingExpression(pheno1, counts1), expected_output)
})

# Test case 2: Check output when all samples have expression data
test_that("function does not remove any samples when all samples have expression data", {
    counts2 <- data.frame(A = c(1, 2, 3, 4, 5), B = c(4, 5, 6, 7, 8), C = c(9, 10, 11, 12, 13))
    colnames(counts2) <- c("sample1", "sample2", "sample3")
    pheno2 <- data.frame(id = c("sample1", "sample2", "sample3"),
                         group = c("control", "treatment", "treatment"))
    rownames(pheno2) <- pheno2$id
    expect_equal(missingExpression(pheno2, counts2), pheno2)
})

# Test case 3: Check output when none of the samples have expression data
test_that("function removes all samples when none of the samples have expression data", {
    counts3 <- data.frame(A = c(1, 2, 3, 4, 5), B = c(4, 5, 6, 7, 8), C = c(9, 10, 11, 12, 13))
    colnames(counts3) <- c("sample4", "sample5", "sample6")
    pheno3 <- data.frame(id = c("sample1", "sample2", "sample3"),
                         group = c("control", "treatment", "treatment"))
    rownames(pheno3) <- pheno3$id
    expect_error(missingExpression(pheno3, counts3))
})
