#' Plot a scaffoldSpace object
#'
#' This function plot a scaffoldSpace object created by \code{buildScaffold} function.
#' @param space A scaffoldSpace object returned by buildScaffold function.
#' @param title Title for the plot
#' @noRd
#' @return ggplot object.


plotScaffold <- function(space,title){
  model <- space@model
  dims <- space@dims

  # calculate variance explained and add percentage to axis labels
  var_sum <- sum(model$sdev^2)
  var1 <- round(model$sdev[dims[1]]^2/var_sum*100,2)
  var2 <- round(model$sdev[dims[2]]^2/var_sum*100,2)
  xlabel <- paste0( "PC",as.character(dims[1]), " (",as.character(var1), "%", ")" )
  ylabel <- paste0( "PC",as.character(dims[2]), " (",as.character(var2), "%", ")" )

  # Prepare a dataframe for ggplot2
  PC1 <- model$x[,dims[1]]
  PC2 <- model$x[,dims[2]]
  Scaffold_group <- space@label
  df <- data.frame(PC1,PC2,Scaffold_group)


  # calculate centroids
  centroids_df <- suppressMessages(df %>% dplyr::group_by(Scaffold_group)
                                   %>% dplyr::summarise(mean_PC1=mean(PC1),mean_PC2=mean(PC2)))

  # define color scheme
  total_types <- length(unique(Scaffold_group))
  if (total_types>12){
    my_col <- RColorBrewer::brewer.pal(12, "Paired")
    pal <- grDevices::colorRampPalette(my_col)
    my_col <- pal(total_types)
    if (space@plot_mode=="dot"){
      message("More than 12 cell types are to be displayed. Setting 'plot_mode='tiny_label' may yield better visualization.")
    }
  } else{
    my_col <- RColorBrewer::brewer.pal(total_types,"Paired")
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
