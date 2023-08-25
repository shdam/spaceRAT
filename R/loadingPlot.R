#' Generate a loading plot of scaffold PCA
#'
#' This function takes a previously built scaffold space
#' and generates a loading plot to visualize the genes
#' that contribute most to the selected principle components.
#'
#' @inheritParams plotScaffold
#' @inheritParams projectSample
#' @param num_genes Optional parameter indicating number of genes to be shown.
#' @param gene_name Character indicating which type of gene name to show by
#' the side of arrows. By default HGNC symbol is shown.
#  Other options are "ensembl_gene", "ensembl_transcript", "entrez",
#  and "refseq_mrna".
#' @param angle A number for the degree of rotation of labels on the plot.
#' @param df_only The data.frame for the loading plot will be returned
#' instead of the plot.
#' @import ggplot2
#' @export
#' @return A data frame indicating the loading scores of the genes
#' that contribute most to the selected principle components.
#' The loading plot is printed automatically, thus not returned.
#' @usage
#' loadingPlot(
#'     space,
#'     dims = c(1, 2),
#'     num_genes = 3,
#'     gene_name = "hgnc_symbol",
#'     angle = 30,
#'     df_only = FALSE
#'     )
#' @examples
#' utils::data("exprs_dmap", "pData_dmap", package = "spaceRATScaffolds")
#' space <- buildScaffold(exprs_dmap, pData_dmap, "cell_types", data = "logged")
#' loadingPlot(space)
loadingPlot <- function(
        space,
        dims = c(1, 2),
        num_genes = 3,
        gene_name = "hgnc_symbol",
        angle = 30,
        df_only = FALSE){

    pca <- space$pca
    pc1 <- dims[1]; pc2 <- dims[2]

    # calculate variance explained and add percentage to axis labels
    var_sum <- sum(pca$sdev^2)
    var1 <- round(pca$sdev[pc1]^2/var_sum*100,2)
    var2 <- round(pca$sdev[pc2]^2/var_sum*100,2)
    xlabel <- paste0(
        "PC", as.character(pc1), " (", as.character(var1), "%", ")" )
    ylabel <- paste0(
        "PC", as.character(pc2), " (", as.character(var2), "%", ")" )
    # determine most important genes
    pc1 <- paste0("PC", pc1); pc2 <- paste0("PC", pc2)
    datapc <- as.data.frame(pca$rotation[,dims])
    datapc <- convertGeneName(datapc,to=gene_name)
    df <- extractGenes(datapc, num_genes)

    xmin <- min(df[[pc1]]); xmax <- max(df[[pc1]])
    ymin <- min(df[[pc2]]); ymax <- max(df[[pc2]])

# ggplot
    if(df_only) return(df)
    g <- ggplot2::ggplot(data=df)+
        ggplot2::aes(color=.data$class) +
        ggplot2::geom_text(ggplot2::aes(
            x = .data[[pc1]], .data[[pc2]],
            label=.data$gene),size = 3, angle=angle,show.legend = FALSE
            )+
        ggplot2::geom_segment(ggplot2::aes(
            x=0, y=0, xend=.data[[pc1]], yend=.data[[pc2]],
            color=.data$class), arrow=ggplot2::arrow(
                length=ggplot2::unit(0.1,"cm")), alpha=0.75
            )+
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

#' Helper function that extracts the relevant number of genes
#' @noRd
extractGenes <- function(datapc, num_genes){
    pc1_ordered <- datapc[order(datapc[,1], decreasing = TRUE),]
    pc2_ordered <- datapc[order(datapc[,2], decreasing = TRUE),]
    pos1 <- data.frame(
        utils::head(pc1_ordered,num_genes),
        class=rep(paste0(
            "Top ",num_genes," genes for ", colnames(datapc)[1]),num_genes))
    neg1 <- data.frame(
        utils::tail(pc1_ordered,num_genes),
        class=rep(paste0(
            "Top ",num_genes," genes for ", colnames(datapc)[1]),num_genes))
    pos2 <- data.frame(
        utils::head(pc2_ordered,num_genes),
        class=rep(paste0(
            "Top ",num_genes," genes for ", colnames(datapc)[2]),num_genes))
    neg2 <- data.frame(
        utils::tail(pc2_ordered,num_genes),
        class=rep(paste0(
            "Top ",num_genes," genes for ", colnames(datapc)[2]),num_genes))
    df <- rbind(pos1,neg1,pos2,neg2)
    df$gene <- rownames(df)
    return(df)
}
