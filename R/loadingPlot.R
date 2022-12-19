#' Generate a loading plot of scaffold PCA
#'
#' This function takes a \code{\link{scaffoldSpace}} object and generates a loading plot to visualize the genes that contribute most to the selected principle components.
#'
#' @param space A scaffoldSpace object created by function \code{\link{buildScaffold}}.
#' @param num_genes An optional parameter indicating number of genes to be shown.
#' @param gene_name A character indicating which type of gene name to show by the side of arrows. By default HGNC symbol is shown.
#  Other options are "ensembl_gene", "ensembl_transcript", "entrez" and "refseq_mrna".
#' @param angle A number for the degree of rotation of labels on the plot.
#' @param df_only The dataframe for the loading plot will be returned instead of the plot.
#'
#' @export
#' @return A data frame indicating the loading scores of the genes that contribute most to the selected principle components.
#' The loading plot is printed automatically, thus not returned.
#'
#' @examples
#' utils::data("exprs_dmap", "pData_dmap", package = "spaceRAT")
#' space <- buildScaffold(exprs_dmap, pData_dmap, "cell_types")
#' loadingPlot(space)

loadingPlot <- function(space,num_genes=3,gene_name="hgnc_symbol",angle=30, df_only=FALSE){
    stopifnot("Please only use loading plot with PCA scaffolds" = is(space@model, "prcomp"))
    pca <- space@model
    pcs <- space@dims

    # calculate variance explained and add percentage to axis labels
    var_sum <- sum(pca$sdev^2)
    var1 <- round(pca$sdev[pcs[1]]^2/var_sum*100,2)
    var2 <- round(pca$sdev[pcs[2]]^2/var_sum*100,2)
    xlabel <- paste0( "PC",as.character(pcs[1]), " (",as.character(var1), "%", ")" )
    ylabel <- paste0( "PC",as.character(pcs[2]), " (",as.character(var2), "%", ")" )

    # determine most important PCs
    datapc <- as.data.frame(pca$rotation[,pcs])
    datapc <- convertGeneName(datapc,to=gene_name)

    pc1_ordered <- datapc[order(datapc[,1], decreasing = T),]
    pc2_ordered <- datapc[order(datapc[,2], decreasing = T),]
    pos1 <- data.frame(utils::head(pc1_ordered,num_genes),class=rep(paste0("Top ",num_genes," genes for PC",as.character(pcs[1])),num_genes))
    neg1 <- data.frame(utils::tail(pc1_ordered,num_genes),class=rep(paste0("Top ",num_genes," genes for PC",as.character(pcs[1])),num_genes))
    pos2 <- data.frame(utils::head(pc2_ordered,num_genes),class=rep(paste0("Top ",num_genes," genes for PC",as.character(pcs[2])),num_genes))
    neg2 <- data.frame(utils::tail(pc2_ordered,num_genes),class=rep(paste0("Top ",num_genes," genes for PC",as.character(pcs[2])),num_genes))
    df <- rbind(pos1,neg1,pos2,neg2) %>% tibble::rownames_to_column(var=gene_name)

    xmin <- min(df[2])
    xmax <- max(df[2])
    ymin <- min(df[3])
    ymax <- max(df[3])

# ggplot
    if(df_only) return(df)
    g <-ggplot2::ggplot(data=df)+
        ggplot2::geom_text(ggplot2::aes_string(colnames(df)[2], colnames(df)[3],label=colnames(df)[1],color=colnames(df)[4]),size = 3, angle=angle,show.legend = FALSE)+
        ggplot2::geom_segment(ggplot2::aes_string(x=0, y=0, xend=colnames(df)[2], yend=colnames(df)[3],color=colnames(df)[4]), arrow=ggplot2::arrow(length=ggplot2::unit(0.1,"cm")), alpha=0.75)+
        ggplot2::theme_bw()+
        ggplot2::theme(legend.text=ggplot2::element_text(hjust=1))+
        ggplot2::xlim(xmin*1.5,xmax*1.5)+
        ggplot2::ylim(ymin*1.5,ymax*1.5)+
        ggplot2::xlab(xlabel)+ggplot2::ylab(ylabel)+
        ggplot2::coord_fixed()+
        ggplot2::ggtitle("Loading plot") +
        ggplot2::labs(color = "Class")


    return(g)#df = df[order(df[,3],decreasing=T),])
}

