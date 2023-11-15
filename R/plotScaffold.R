#' Plot a scaffoldSpace object
#'
#' This function plot a \code{scaffoldSpace} object created by
#' \code{\link{buildScaffold}} function.
#' @param scaffold A \code{scaffoldSpace} object returned by
#' \code{\link{buildScaffold}} function.
#' @param title Title for the plot
#' @param plot_mode A character indicating whether to add tiny
#' By default \code{plot_mode="dot"} and tiny labels will not be attached.
#' If more than 12 cell types are to be displayed, setting
#' \code{plot_mode="tiny_label"} may yield better visualization.
#' Shorter names for phenotypes (e.g. cell types) is strongly recommended
#' in "tiny_label" mode.
#' @param dimred A character indicating the method for
#' dimensionality reduction. Currently "PCA" and "UMAP" are supported.
#' labels to each data point.
#' @param dims A numeric vector containing 2 numbers, indicating
#' which two principle components to plot.
#' @export
#' @importFrom stats aggregate
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @import ggplot2
#' @return ggplot object.
#' @usage
#' plotScaffold(
#' scaffold,
#' dimred = "PCA",
#' title = "Scaffold plot",
#' plot_mode = "dot",
#' dims = c(1, 2))
#' @examples
#' scaffold <- buildScaffold("DMAP_scaffold")
#' plotScaffold(scaffold, "Scaffold plot title", dimred = "PCA")
plotScaffold <- function(
        scaffold,
        dimred = "PCA",
        title = "Scaffold plot",
        plot_mode = "dot",
        dims = c(1, 2)){

    stopifnot("Please use either 'PCA' or 'UMAP' as dimred." = toupper(dimred) %in% c("PCA", "UMAP"))

    if(toupper(dimred) == "PCA"){
        stopifnot("No PCA space in scaffold" = !is(scaffold$pca, "NULL"))
        # calculate variance explained and add percentage to axis labels
        var_sum <- sum(scaffold$pca$sdev^2)
        var1 <- round(scaffold$pca$sdev[dims[1]]^2/var_sum*100,2)
        var2 <- round(scaffold$pca$sdev[dims[2]]^2/var_sum*100,2)
        xlabel <- paste0(
            "PC",as.character(dims[1]),
            " (",as.character(var1), "%", ")" )
        ylabel <- paste0(
            "PC",as.character(dims[2]),
            " (",as.character(var2), "%", ")" )

        # Prepare a dataframe for ggplot2
        Dim1 <- scaffold$pca$x[,dims[1]]
        Dim2 <- scaffold$pca$x[,dims[2]]
    } else if(toupper(dimred) == "UMAP"){
        stopifnot("No UMAP space in scaffold" = !is(scaffold$umap, "NULL"))
        # Prepare a dataframe for ggplot2
        Dim1 <- scaffold$umap$embedding[,1]
        Dim2 <- scaffold$umap$embedding[,2]
        xlabel <- "UMAP1"
        ylabel <- "UMAP2"
    }


    df <- data.frame(Dim1,Dim2,"Scaffold_group" = scaffold$label)


    # calculate centroids
    centroids_df <- stats::aggregate(
        cbind(Dim1, Dim2) ~ df$Scaffold_group, FUN = mean)
    colnames(centroids_df) <- c("Scaffold_group", "mean_Dim1", "mean_Dim2")


    # define color scheme
    total_types <- length(unique(df$Scaffold_group))
    if (total_types>12){
        my_col <- RColorBrewer::brewer.pal(12, "Paired")
        pal <- grDevices::colorRampPalette(my_col)
        my_col <- pal(total_types)
        if (plot_mode=="dot"){
            message(
            "More than 12 cell types are to be displayed. Setting
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
                x = .data$mean_Dim1,
                y = .data$mean_Dim2,
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
        message(
        "plot_mode='tiny_label', shorter names for cell types in
        phenotype table yields better visualization.")
        g <- ggplot2::ggplot()+
            ggplot2::geom_point(data=df,mapping=ggplot2::aes(
                x = .data$Dim1,
                y = .data$Dim2,
                color=.data$Scaffold_group),
                size = 0.8)+
            ggplot2::geom_text(
                data=df,ggplot2::aes(
                    x = .data$Dim1, y = .data$Dim2,
                    label = .data$Scaffold_group, color = .data$Scaffold_group),
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
