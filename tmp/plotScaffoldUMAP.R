#' Plot a scaffoldSpace object containing a UMAP model.
#'
#' This function plot a scaffoldSpace object containing a UMAP model.
#' @param space A scaffoldSpace object containing a UMAP model.
#' @param title Title for the plot.
#' @noRd
#' @return a ggplot object.


plotScaffoldUMAP<- function(space,title){

    # Prepare a dataframe for ggplot2
    Dim1 <- space@model[,1]
    Dim2 <- space@model[,2]
    Scaffold_group <- space@label
    df <- data.frame(Dim1,Dim2,Scaffold_group)


    # calculate centroids
    centroids_df <- suppressMessages(
        df %>%
            dplyr::group_by(Scaffold_group) %>%
            dplyr::summarise(mean_Dim1 = mean(Dim1),
                             mean_Dim2 = mean(Dim2)))

    # define color scheme
    total_types <- length(unique(Scaffold_group))
    my_col <- RColorBrewer::brewer.pal(total_types,"Paired")
    if (total_types>12){
        pal <- grDevices::colorRampPalette(my_col)
        my_col <- pal(total_types)
        if (space@plot_mode=="dot"){
          warning("More than 12 cell types are to be displayed. Setting 'plot_mode='tiny_label' may yield better visualization.")
        }
    }

    # ggplot2 for dot mode
    if (space@plot_mode=="dot"){
        suppressMessages(
            g <- ggplot2::ggplot(data=df)+
                ggplot2::geom_point(mapping=ggplot2::aes(Dim1,Dim2,color=Scaffold_group))+
                ggplot2::scale_color_manual(name="Scaffold_group",values=my_col)+
                ggplot2::geom_label(data=centroids_df,ggplot2::aes(mean_Dim1,mean_Dim2,label=Scaffold_group,color=Scaffold_group),fontface="bold",show.legend = FALSE)+
                ggplot2::labs(
                    title = title,
                    x = "UMAP1",
                    y = "UMAP2"
                ) +
                ggplot2::theme_bw()+
                ggplot2::coord_fixed()
        )
        return(g)
    }

    # ggplot2 for label mode
    if (space@plot_mode=="tiny_label"){
        message("plot_mode='tiny_label', shorter names for cell types in phenotype table yields better visualization.")
        g <- ggplot2::ggplot()+
            ggplot2::geom_point(data=df,mapping=ggplot2::aes(Dim1,Dim2,color=Scaffold_group),size = 0.8)+
            ggplot2::geom_text(data=df,ggplot2::aes(Dim1,Dim2,label=Scaffold_group,color=Scaffold_group),alpha=0.5,size=3,show.legend = FALSE)+
            ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=2)))+
            ggplot2::labs(
                title = title,
                x = "UMAP1",
                y = "UMAP2"
            ) +
            ggplot2::scale_color_manual(name="Scaffold_group",values=my_col)+
            ggplot2::theme_bw()+
            ggplot2::coord_fixed()
        return(g)
    }
    }
