#' Plot a scaffoldSpace object
#'
#' This function plot a scafoldSpace object created by \code{buildScaffold} function.
#' @import magrittr
#' @import  ggplot2
#' @importFrom dplyr filter group_by summarise
#' @importFrom Biobase pData exprs
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @param space A scaffoldSpace object returned by buildScaffold function.
#' @param title Title for the plot
#' @param classes Classes not to show.
#' @noRd
#' @return NULL. This function doesn't return any value.


plotScaffold <- function(space,title,classes=NULL){
        pca <- space@pca
        pcs <- space@pcs

        # calculate variance explained and add percentage to axis labels
        var_sum <- sum(pca$sdev^2)
        var1 <- round(pca$sdev[pcs[1]]^2/var_sum*100,2)
        var2 <- round(pca$sdev[pcs[2]]^2/var_sum*100,2)
        xlabel <- paste0( "PC",as.character(pcs[1]), " (",as.character(var1), "%", ")" )
        ylabel <- paste0( "PC",as.character(pcs[2]), " (",as.character(var2), "%", ")" )

        # Prepare a dataframe for ggplot2
        PC1 <- pca$x[,pcs[1]]
        PC2 <- pca$x[,pcs[2]]
        Scaffold_group <- space@label
        df <- data.frame(PC1,PC2,Scaffold_group)

        if(!is.null(classes)) {
                df <- df %>%
                      dplyr::filter(class %in% classes)
        }

        # calculate centroids
        centroids_df <- suppressMessages(df %>% dplyr::group_by(Scaffold_group)
                                         %>% dplyr::summarise(mean_PC1=mean(PC1),mean_PC2=mean(PC2)))

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
                     ggplot2::geom_point(mapping=ggplot2::aes(PC1,PC2,color=Scaffold_group))+
                     ggplot2::scale_color_manual(name="Scaffold_group",values=my_col)+
                     ggplot2::geom_label(data=centroids_df,ggplot2::aes(mean_PC1,mean_PC2,label=Scaffold_group,color=Scaffold_group),fontface="bold",show.legend = FALSE)+
                     ggplot2::ggtitle(title)+
                     ggplot2::xlab(xlabel)+
                     ggplot2::ylab(ylabel)+
                     ggplot2::theme_bw()+
                     ggplot2::coord_fixed()
                )
                return(g)
        }

  # ggplot2 for label mode
        if (space@plot_mode=="tiny_label"){
                message("plot_mode='tiny_label', shorter names for cell types in phenotype table yields better visualization.")
                g <- ggplot2::ggplot()+
                     ggplot2::geom_point(data=df,mapping=ggplot2::aes(PC1,PC2,color=Scaffold_group),size = 0.8)+
                     ggplot2::geom_text(data=df,ggplot2::aes(PC1,PC2,label=Scaffold_group,color=Scaffold_group),alpha=0.5,size=3,show.legend = FALSE)+
                     ggplot2::geom_label(data=centroids_df,ggplot2::aes(mean_PC1,mean_PC2,label=Scaffold_group,color=Scaffold_group,fontface="bold"),show.legend = FALSE)+
                     ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=2)))+
                     ggplot2::ggtitle(title)+
                     ggplot2::xlab(xlabel)+
                     ggplot2::ylab(ylabel)+
                     ggplot2::scale_color_manual(name="Scaffold_group",values=my_col)+
                     ggplot2::theme_bw()+
                     ggplot2::coord_fixed()
                return(g)
        }

}
