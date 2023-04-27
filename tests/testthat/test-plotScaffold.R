test_that("plotScaffold works", {
    scaffold <- buildScaffold("prebuilt_DMAP")
    plt <- plotScaffold(scaffold, "Scaffold plot title")
    expect_s3_class(plt, "ggplot")

    scaffold@plot_mode <- "tiny_label"
    plt <- plotScaffold(scaffold, "Scaffold plot title")
    expect_s3_class(plt, "ggplot")

    scaffold@model <- scaffold@model$x[,1:2]
    scaffold@plot_mode <- "dot"
    plt <- plotScaffold(scaffold, "Scaffold plot title")
    expect_s3_class(plt, "ggplot")

    scaffold@label <- paste0(scaffold@label, 1:2)
    scaffold@plot_mode <- "tiny_label"
    plt <- plotScaffold(scaffold, "Scaffold plot title")
    expect_s3_class(plt, "ggplot")

    scaffold@plot_mode <- "dot"
    plt <- plotScaffold(scaffold, "Scaffold plot title")
    expect_s3_class(plt, "ggplot")
})
