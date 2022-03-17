test_that("projectSample() correctly projects.",{
  space <- buildScaffold("prebuilt_DMAP")
  g <-projectSample(space,exprs_ilaria,pData_ilaria,"cancer_type")
  expect_s3_class(g,"ggplot")
})
