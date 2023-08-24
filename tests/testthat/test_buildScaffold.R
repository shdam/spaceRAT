data("exprs_dmap", package = "spaceRATScaffolds")
data("pData_dmap", package = "spaceRATScaffolds")

test_that("buildScaffold() returns an scaffoldSpace object",{

    scaffold <- "DMAP_scaffold"
    # test prebuilt scaffold
    space1 <- buildScaffold(scaffold, data = "logged")
    expect_type(space1,"list")
    expect_setequal(names(space1), c("label","DEgenes", "pca"))

    space2 <- buildScaffold(scaffold,classes=c("HSC","MONO","ERY"), data = "logged")
    expect_type(space2,"list")
    expect_setequal(names(space2), c("label","DEgenes", "pca"))
    expect_setequal(as.character(unique(space2$label)), c("HSC","MONO","ERY"))

    space3 <- buildScaffold(scaffold,classes=c("BASO","EOS","GRAN"), data = "logged")
    expect_type(space3,"list")
    expect_setequal(names(space3), c("label","DEgenes", "pca"))
    expect_setequal(as.character(unique(space3$label)), c("BASO","EOS","GRAN"))

    space4 <- buildScaffold(scaffold,classes=c("MEGA","ERY","MEP","HSC"), data = "logged")
    expect_type(space4,"list")
    expect_setequal(names(space4), c("label","DEgenes", "pca"))
    expect_setequal(as.character(unique(space4$label)), c("MEGA","ERY","MEP","HSC"))

    # test newly-built scaffold
    space5 <- buildScaffold(exprs_dmap, pData_dmap, "cell_types", data = "logged")
    expect_type(space5,"list")
    expect_setequal(names(space5), c("label","DEgenes", "pca"))

    space6 <- buildScaffold(exprs_dmap,pData_dmap,"cell_types",classes=c("HSC","MEP","ERY"), data = "logged")
    expect_type(space6,"list")
    expect_setequal(names(space6), c("label","DEgenes", "pca"))
    expect_setequal(as.character(unique(space6$label)), c("HSC","MEP","ERY"))

    space7 <- buildScaffold(exprs_dmap, pData_dmap, data = "logged", add_umap = TRUE)
    expect_type(space7,"list")
    expect_setequal(names(space7), c("label","DEgenes", "pca", "umap"))

})

test_that("test error checks", {

    expect_error(buildScaffold("something"), "Incorrectly named prebuilt scaffold. The available are.*")

    colnames(exprs_dmap) <- pData_dmap$cell_types
    expect_error(buildScaffold(exprs_dmap), "Please ensure unique column names in data.")


})
