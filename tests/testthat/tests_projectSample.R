data("counts_ilaria", package = "spaceRATScaffolds")
data("pData_ilaria", package = "spaceRATScaffolds")
space <- buildScaffold("DMAP_scaffold")
test_that("projectSample() correctly projects.",{
    g <- projectSample(space,counts_ilaria,pData_ilaria,"cancer_type")
    expect_s3_class(g,"ggplot")

    g <- projectSample(space,counts_ilaria)
    expect_s3_class(g,"ggplot")
})

test_that("projectSample() gives warning when sample is shallow.",{
    expect_warning(projectSample(space,counts_ilaria[1:1000,1:10]), "More than 1/4 genes are added to sample with imputed.*")
})
