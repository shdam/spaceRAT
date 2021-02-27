#' Returns a ggplot and PCA parameters for ExpressionSet object
#'
#' This function performs differential expression analysis to select DE gene, then subset the scaffold dataset, rank the subset, finally produce a PCA plot using ranks.
#' @import ggplot2
#' @importFrom Biobase pData
#' @param exprs_scaffold an expression matrix of class matrix, or data frame that can be converted to matrix.
#' @param pData_scaffold a dataframe corresponding to the expression matrix. The row names of pData_scaffold must be identical to column names of exprs_scaffold.
#' @param group_scaffold a column name of pData_scaffold
#' @param pval_cutoff a cutoff value for p value when selecting differentially expressed genes. By default 0.05
#' @param lfc_cutoff a cutoff value for logFC when selecting differentially expressed genes. By default 2
#' @param pca_scale a logical variable determining whether to normalize rows when plotting PCA
#' @export
#' @return A scaffoldSpace object
#' @examples
#' buildScaffold(exprs_dmap,pData_dmap"cell_types")
#' buildScaffold(exprs_dmap,pData_dmap,"cell_types", pval_cutoff=0.01,pca_scale=T)

buildScaffold <- function(exprs_scaffold,
                          pData_scaffold,
                          group_scaffold,
                          classes = NULL,
                          pval_cutoff=0.05,
                          lfc_cutoff=2,
                          title="Scaffold PCA Plot",
                          pca_scale=FALSE){
        # create eset
        eset_scaffold <- createEset(exprs_scaffold,pData_scaffold)

        # find DE genes
        DEgenes <- findDEGenes(eset_scaffold,group_scaffold,pval_cutoff,lfc_cutoff)

        # subset
        eset_scaffold <- eset_scaffold[DEgenes,]

        # rank
        scaffold_rank<- apply(Biobase::exprs(eset_scaffold),2,rank)

        # PCA
        pca <- prcomp(t(scaffold_rank),scale=pca_scale)

        # Prepare a dataframe for ggplot2
        PC1 <- pca$x[,1]
        PC2 <- pca$x[,2]
        scaffold_group <- Biobase::pData(eset_scaffold)[[group_scaffold]]
        df <- data.frame(PC1,PC2,scaffold_group)
        if(!is.null(classes)) {
          df <- df %>%
            dplyr::filter(scaffold_group %in% classes)
        }

        # ggplot2
        g <- ggplot2::ggplot()+ggplot2::geom_point(data=df,mapping=ggplot2::aes(PC1,PC2,color=scaffold_group))+ggplot2::ggtitle(title)
        # print(g)

        # record standard data in scaffold
        scaffold <- new("scaffoldSpace",
                        graph=g,
                        pca=pca,
                        DEgenes=DEgenes)

        return(scaffold)
}
