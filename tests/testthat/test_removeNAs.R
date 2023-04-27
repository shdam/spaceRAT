test_that("function replaces NAs with 0", {
    counts1 <- c(1, 2, NA, 4)
    counts2 <- c(NA, 3, 4, 5)
    counts3 <- c(1, 2, 3, 4)

    expect_equal(removeNAs(counts1), c(1, 2, 0, 4))
    expect_equal(removeNAs(counts2), c(0, 3, 4, 5))
    expect_equal(removeNAs(counts3), c(1, 2, 3, 4))
})

test_that("function gives proper message when there are NAs", {
    counts1 <- c(1, 2, NA, 4)
    counts2 <- c(NA, 3, 4, 5)

    expect_message(removeNAs(counts1),
                   "Expression data has 1 missing values. Replacing NAs by 0.")
    expect_message(removeNAs(counts2),
                   "Expression data has 1 missing values. Replacing NAs by 0.")
})

test_that("function gives no message when there are no NAs", {
    counts1 <- c(1, 2, 3, 4)

    expect_equal(removeNAs(counts1), counts1)
})
