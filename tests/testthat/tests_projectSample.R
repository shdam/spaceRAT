data("ilaria_counts", package = "spaceRATScaffolds")
data("ilaria_pData", package = "spaceRATScaffolds")
space <- buildScaffold("DMAP.v1")
test_that("projectSample() correctly projects.",{
    g <- projectSample(space,ilaria_counts,ilaria_pData,"cancer_type")
    expect_s3_class(g,"ggplot")

    g <- projectSample(space,ilaria_counts,ilaria_pData,"cancer_type", dims = c(2,3))
    expect_s3_class(g,"ggplot")

    g <- projectSample(space,ilaria_counts,ilaria_pData,"cancer_type", title = "Some title")
    expect_s3_class(g,"ggplot")

    g <- projectSample(space,ilaria_counts)
    expect_s3_class(g,"ggplot")

    g <- projectSample(space,ilaria_counts, dims = c(2,3))
    expect_s3_class(g,"ggplot")

})

test_that("projectSample() gives warning when sample is shallow.",{
    expect_warning(projectSample(space,ilaria_counts[1:1000,1:10]), "More than 1/4 genes are added to sample with imputed.*")
})
