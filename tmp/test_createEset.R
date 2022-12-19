test_that("createEset() returns an ExpressionSet object with proper dimension",{
  eset_dmap <- createEset(exprs_dmap[1:20,1:10],pData_dmap[1:10,,drop=FALSE],"cell_types")

  expect_s4_class(eset_dmap,"ExpressionSet")
  expect_equal(unname(dim(eset_dmap)[1]),20)
  expect_equal(unname(dim(eset_dmap)[2]),10)
})


test_that("createEset() handles NA in expression matrix or phenotype table",{
  exprs_dmap_na <- exprs_dmap[1:20,1:10]
  exprs_dmap_na[1,] <- NA
  exprs_dmap_na[5,5:10] <- NA
  pData_dmap_na <- pData_dmap[1:10,,drop=FALSE]
  pData_dmap_na[2:4,1] <-NA
  res <- evaluate_promise(createEset(exprs_dmap_na,pData_dmap_na,"cell_types"))

  expect_equal(res$messages[1],"3 row(s) in phenotype table contain NA in the required column, thus removed.\n")
  expect_equal(res$messages[2],"3 sample(s) have no phenotype annotation, thus removed.\n")
  expect_equal(res$messages[3],"Count matrix has 13 missing values. Replace NA by 0.\n")

  expect_s4_class(res$result,"ExpressionSet")
  expect_equal(unname(dim(res$result)[1]),20)
  expect_equal(unname(dim(res$result)[2]),7)
})


test_that("createEset() calculates intersection of columns of expression matrix and rows of phenotype table, then make them match.",{
  exprs_messy <- exprs_dmap[1:20,sample(1:15)]
  pData_messy <- pData_dmap[sample(1:10),,drop=FALSE]
  res <- evaluate_promise(createEset(exprs_messy,pData_messy,"cell_types"))
  expect_s4_class(res$result,"ExpressionSet")
  expect_equal(unname(dim(res$result)[1]),20)
  expect_equal(unname(dim(res$result)[2]),10)

  exprs_messy <- exprs_dmap[1:20,sample(1:15)]
  pData_messy <- pData_dmap[sample(1:30),,drop=FALSE]
  res <- evaluate_promise(createEset(exprs_messy,pData_messy,"cell_types"))
  expect_s4_class(res$result,"ExpressionSet")
  expect_equal(unname(dim(res$result)[1]),20)
  expect_equal(unname(dim(res$result)[2]),15)
})


test_that("createEset() properly subset",{
  eset <- createEset(exprs_dmap,pData_dmap, "cell_types",classes =  c("HSC","CMP"))
  expect_equal(unname(dim(eset)[2]),18)

  eset <- createEset(exprs_dmap,pData_dmap,"cell_types",c("ERY","BASO","MONO"))
  expect_equal(unname(dim(eset)[2]),48)
})
