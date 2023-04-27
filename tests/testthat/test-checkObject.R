
# Generate example data
set.seed(123)
counts <- matrix(rpois(24, 10), ncol = 4, dimnames = list(paste0("Gene", 1:6), paste0("Sample", 1:4)))
colData <- data.frame(condition = c(rep("control", 2), rep("treatment", 2)), stringsAsFactors = FALSE)
rownames(colData) <- colnames(counts)

test_that("checkObject returns an object of class matrix or data.frame", {
    obj <- checkObject(counts)
    expect_true(is(obj, "matrix"))
    obj <- checkObject(data.frame(counts))
    expect_s3_class(obj, "data.frame")
})

test_that("checkObject returns a list when input object is not a matrix or data.frame", {
    se <- SummarizedExperiment::SummarizedExperiment(assays = S4Vectors::SimpleList(counts = counts), colData = colData)
    obj <- checkObject(se, assay = "counts")
    obj2 <- checkObject(se, assay = NULL)
    expect_type(obj, "list")
    expect_length(obj, 2)
    expect_identical(obj[[1]], counts)
    expect_identical(obj[[2]], colData)
    expect_identical(obj, obj2)
})

test_that("checkObject fails if wrong format is provided", {
    expect_error(checkObject(list(counts, colData)), "Expression data was not provided in a supported format.*")
    expect_error(checkObject(NULL), "Expression data was not provided in a supported format.*")
})
