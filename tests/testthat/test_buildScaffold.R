test_that("buildScaffold() returns an scaffoldSpace object",{

    data("exprs_dmap")
    data("pData_dmap")
    # test prebuilt scaffold
    space1 <- buildScaffold("prebuilt_DMAP")
    expect_s4_class(space1,"scaffoldSpace")

    space2 <- buildScaffold("prebuilt_DMAP",classes=c("HSC","MONO","ERY"))
    expect_s4_class(space2,"scaffoldSpace")

    space3 <- buildScaffold("prebuilt_DMAP",classes=c("BASO","EOS","GRAN"))
    expect_s4_class(space3,"scaffoldSpace")

    space4 <- buildScaffold("prebuilt_DMAP",classes=c("MEGA","ERY","MEP","HSC"))
    expect_s4_class(space4,"scaffoldSpace")

    # test newly-built scaffold
    space5 <- buildScaffold(exprs_dmap,pData_dmap,"cell_types",plot_mode="tiny_label",dims=c(1,2))
    expect_s4_class(space5,"scaffoldSpace")

    space6 <- buildScaffold(exprs_dmap,pData_dmap,"cell_types",plot_mode="tiny_label",classes=c("HSC","MEP","ERY"),dims=c(1,2))
    expect_s4_class(space6,"scaffoldSpace")

    space7 <- buildScaffold(exprs_dmap, pData_dmap, dim_reduction = "UMAP")
    expect_s4_class(space7,"scaffoldSpace")

})

test_that("test warnings", {

    expect_error(buildScaffold("something"), "Incorrectly named prebuilt scaffold. The available are.*")

    colnames(exprs_dmap) <- pData_dmap$cell_types
    expect_warning(buildScaffold(exprs_dmap), "No annotation data provided.*")


})
