data("exprs_ilaria")
data("pData_ilaria")
space <- buildScaffold("prebuilt_DMAP")
test_that("projectSample() correctly projects.",{
    g <- projectSample(space,exprs_ilaria,pData_ilaria,"cancer_type")
    expect_s3_class(g,"ggplot")

    g <- projectSample(space,exprs_ilaria)
    expect_s3_class(g,"ggplot")
})

test_that("projectSample() gives warning when sample is shallow.",{
    expect_warning(projectSample(space,exprs_ilaria[1:1000,1:10]), "More than 1/4 genes are added to sample with imputed.*")
})
