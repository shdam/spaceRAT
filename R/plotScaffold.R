#' Plot a scaffoldSpace object
#'
#' This function plot a \code{scaffoldSpace} object created by
#' \code{\link{buildScaffold}} function.
#' @param space A \code{scaffoldSpace} object returned by
#' \code{\link{buildScaffold}} function.
#' @param title Title for the plot
#' @param plot_mode A character indicating whether to add tiny
#' labels to each data point.
#' By default \code{plot_mode="dot"} and tiny labels will not be attached.
#' If more than 12 cell types are to be displayed, setting
#' \code{plot_mode="tiny_label"} may yield better visualization.
#' Shorter names for phenotypes (e.g. cell types) is strongly recommended
#' in "tiny_label" mode.
#' @param dims A numeric vector containing 2 numbers, indicating
#' which two principle components to plot.
#' @export
#' @importFrom stats aggregate
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @import ggplot2
#' @return ggplot object.
#' @examples
#' scaffold <- buildScaffold("DMAP_scaffold")
#' plotScaffold(scaffold, "Scaffold plot title")
plotScaffold <- function(
        space,
        dim_reduction = c("PCA", "UMAP"),
        title = "Scaffold plot",
        plot_mode = c("dot", "tiny_label"),
        dims = c(1, 2)){

    plot_mode <- match.arg(plot_mode)
    dim_reduction <- match.arg(dim_reduction)

    if(dim_reduction == "PCA"){
        stopifnot("No PCA space in scaffold" = !is(space$pca, "NULL"))
        # calculate variance explained and add percentage to axis labels
        var_sum <- sum(space$pca$sdev^2)
        var1 <- round(space$pca$sdev[dims[1]]^2/var_sum*100,2)
        var2 <- round(space$pca$sdev[dims[2]]^2/var_sum*100,2)
        xlabel <- paste0(
            "PC",as.character(dims[1]),
            " (",as.character(var1), "%", ")" )
        ylabel <- paste0(
            "PC",as.character(dims[2]),
            " (",as.character(var2), "%", ")" )

        # Prepare a dataframe for ggplot2
        Dim1 <- space$pca$x[,dims[1]]
        Dim2 <- space$pca$x[,dims[2]]
    } else{
        stopifnot("No UMAP space in scaffold" = !is(space$umap, "NULL"))
        # Prepare a dataframe for ggplot2
        Dim1 <- space$umap[,1]
        Dim2 <- space$umap[,2]
        xlabel <- "UMAP1"
        ylabel <- "UMAP2"
    }


    df <- data.frame(Dim1,Dim2,"Scaffold_group" = space$label)


    # calculate centroids
    centroids_df <- stats::aggregate(
        cbind(Dim1, Dim2) ~ Scaffold_group, data = df, FUN = mean)
    colnames(centroids_df) <- c("Scaffold_group", "mean_Dim1", "mean_Dim2")


    # define color scheme
    total_types <- length(unique(df$Scaffold_group))
    if (total_types>12){
        my_col <- RColorBrewer::brewer.pal(12, "Paired")
        pal <- grDevices::colorRampPalette(my_col)
        my_col <- pal(total_types)
        if (plot_mode=="dot"){
            message("More than 12 cell types are to be displayed. Setting
                    'plot_mode='tiny_label' may yield better visualization.")
        }
    } else{
        my_col <- RColorBrewer::brewer.pal(total_types, "Paired")
    }

    # ggplot2 for dot mode
    if (plot_mode=="dot"){
        g <- ggplot2::ggplot(data=df) +
            ggplot2::geom_point(
                mapping=ggplot2::aes(
                    .data$Dim1,.data$Dim2, color=.data$Scaffold_group)) +
            ggplot2::scale_color_manual(name="Scaffold_group",values=my_col) +
            ggplot2::geom_label(data=centroids_df,ggplot2::aes(
                .data$mean_Dim1,
                .data$mean_Dim2,
                label=.data$Scaffold_group,color=.data$Scaffold_group),
                fontface="bold",show.legend = FALSE) +
            ggplot2::labs(
                title = title,
                x = xlabel,
                y = ylabel
            ) +
            ggplot2::theme_bw() +
            ggplot2::coord_fixed()
        return(g)
    }

    # ggplot2 for label mode
    if (plot_mode=="tiny_label"){
        message("plot_mode='tiny_label', shorter names for cell types in
                phenotype table yields better visualization.")
        g <- ggplot2::ggplot()+
            ggplot2::geom_point(data=df,mapping=ggplot2::aes(
                Dim1,Dim2,color=Scaffold_group),size = 0.8)+
            ggplot2::geom_text(data=df,ggplot2::aes(
                Dim1,Dim2,label=Scaffold_group,color=Scaffold_group),
                alpha=0.5,size=3,show.legend = FALSE)+
            ggplot2::geom_label(data=centroids_df,ggplot2::aes(
                .data$mean_Dim1,.data$mean_Dim2,
                label=.data$Scaffold_group,color=.data$Scaffold_group,
                fontface="bold"),show.legend = FALSE)+
            ggplot2::guides(colour = ggplot2::guide_legend(
                override.aes = list(size=2)))+
            ggplot2::labs(
                title = title,
                x = xlabel,
                y = ylabel
            ) +
            ggplot2::scale_color_manual(name="Scaffold_group",values=my_col)+
            ggplot2::theme_bw()+
            ggplot2::coord_fixed()
        return(g)
    }

}
