data("DMAP_exprs", package = "spaceRATScaffolds")
data("DMAP_pData", package = "spaceRATScaffolds")

test_that("buildScaffold() returns a proper space",{

    scaffold <- "DMAP.v1"
    # test prebuilt scaffold
    space1 <- buildScaffold(scaffold, path = "../scaffolds")
    expect_type(space1,"list")
    expect_setequal(names(space1), c("label","DEgenes", "pca", "umap"))

    space2 <- buildScaffold(scaffold,classes=c("HSC","MONO","ERY"), data = "exprs", add_umap = TRUE)
    expect_type(space2,"list")
    expect_true("umap" %in% names(space2))
    expect_true("pca" %in% names(space2))
    expect_setequal(as.character(unique(space2$label)), c("HSC","MONO","ERY"))

    space3 <- buildScaffold(scaffold,classes=c("BASO","EOS","GRAN"), data = "exprs")
    expect_type(space3,"list")
    expect_false("umap" %in% names(space3))
    expect_true("pca" %in% names(space3))
    expect_setequal(as.character(unique(space3$label)), c("BASO","EOS","GRAN"))

    space4 <- buildScaffold(scaffold,classes=c("MEGA","ERY","MEP","HSC"), data = "exprs")
    expect_type(space4,"list")
    expect_false("umap" %in% names(space4))
    expect_true("pca" %in% names(space4))
    expect_setequal(as.character(unique(space4$label)), c("MEGA","ERY","MEP","HSC"))

    # test newly-built scaffold
    space5 <- buildScaffold(DMAP_exprs, DMAP_pData, "cell_types", data = "exprs")
    expect_type(space5,"list")
    expect_true("pca" %in% names(space5))

    space6 <- buildScaffold(DMAP_exprs,DMAP_pData,"cell_types",classes=c("HSC","MEP","ERY"), data = "exprs")
    expect_type(space6,"list")
    expect_true("pca" %in% names(space6))
    expect_setequal(as.character(unique(space6$label)), c("HSC","MEP","ERY"))

    space7 <- buildScaffold(DMAP_exprs, DMAP_pData,"cell_types", data = "exprs", add_umap = TRUE)
    expect_type(space7,"list")
    expect_true("umap" %in% names(space7))
    expect_true("pca" %in% names(space7))

})

test_that("test error checks", {

    expect_error(buildScaffold("something"), "something is not an available scaffold. The available are:*")

    expect_error(buildScaffold(DMAP_exprs, pheno = DMAP_pData, data = "exprs"), "Please specify colname for pheno data")
    expect_error(buildScaffold(DMAP_exprs, pheno = NULL, data = "exprs"), "All cells have unique phenotype information.*")


    colnames(DMAP_exprs) <- DMAP_pData$cell_types
    expect_error(buildScaffold(DMAP_exprs, data = "exprs"), "Please ensure unique column names in data.")


})
