test_that("plotScaffold works", {
    scaffold <- buildScaffold("DMAP_scaffold", classes = c("BASO", "CMP", "EOS", "ERY", "GMP", "GRAN", "HSC", "MEGA", "MEP", "MONO"),
                              add_umap = TRUE, data = "exprs")

    plt <- plotScaffold(scaffold, dimred = "PCA", "Scaffold plot title")
    expect_s3_class(plt, "ggplot")

    plt <- plotScaffold(scaffold, dimred = "UMAP", plot_mode = "tiny_label", "Scaffold plot title")
    expect_s3_class(plt, "ggplot")

    scaffold$pca$x <- scaffold$pca$x[,1:2]
    plt <- plotScaffold(scaffold, dimred = "PCA", plot_mode = "dot", "Scaffold plot title")
    expect_s3_class(plt, "ggplot")

    scaffold$label <- paste0(scaffold$label, 1:2)
    plt <- plotScaffold(scaffold, dimred = "PCA", plot_mode = "tiny_label", "Scaffold plot title")
    expect_s3_class(plt, "ggplot")

    plt <- plotScaffold(scaffold, dimred = "UMAP", plot_mode = "dot", "Scaffold plot title")
    expect_s3_class(plt, "ggplot")
})
