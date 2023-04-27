test_that("function matches row/col names correctly", {
    counts1 <- data.frame(A = c(1, 2, 3), B = c(4, 5, 6), C = c(7, 8, 9))
    colnames(counts1) <- c("sample1", "sample2", "sample3")
    pheno1 <- data.frame(id = c("sample1", "sample2", "sample3"),
                         group = c("control", "treatment", "treatment"))
    rownames(pheno1) <- pheno1$id

    counts2 <- data.frame(A = c(1, 2, 3), B = c(4, 5, 6), C = c(7, 8, 9))
    colnames(counts2) <- c("sample3", "sample4", "sample5")
    pheno2 <- data.frame(id = c("sample1", "sample2", "sample3", "sample4",
                                "sample5"),
                         group = c("control", "treatment", "treatment",
                                   "control", "treatment"))
    rownames(pheno2) <- pheno2$id

    expect_equal(matchToExpression(pheno1, counts1),
                 data.frame(id = c("sample1", "sample2", "sample3"),
                            group = c("control", "treatment", "treatment"),
                            row.names = c("sample1", "sample2", "sample3")))
    expect_equal(matchToExpression(pheno2, counts2),
                 data.frame(id = c("sample3", "sample4", "sample5"),
                 group = c("treatment", "control", "treatment"),
                 row.names = c("sample3", "sample4", "sample5")))
})

test_that("function gives correct output when there is no match", {
    counts1 <- data.frame(A = c(1, 2, 3), B = c(4, 5, 6), C = c(7, 8, 9))
    colnames(counts1) <- c("sample4", "sample5", "sample6")
    pheno1 <- data.frame(id = c("sample1", "sample2", "sample3"),
                         group = c("control", "treatment", "treatment"))
    rownames(pheno1) <- pheno1$id

    expect_equal(matchToExpression(pheno1, counts1),
                 data.frame(id = character(0),
                            group = character(0),
                            row.names = character(0)))
})
