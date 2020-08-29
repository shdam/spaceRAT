#' @title Returns a ggplot and PCA parameters for ExpressionSet object
#' @description This function performs differential expression analysis to select DE gene, then subset the scaffold dataset, rank the subset, finally produce a PCA plot using ranks.
#' @import ggplot2
#' @import Biobase
#' @param eset_scaffold an ExpressionSet
#' @param group_scaffold a column name of pData(eset_scaffold)
#' @param pval_cutoff a cutoff value for p value when selecting differentially expressed genes. By default 0.05
#' @param lfc_cutoff a cutoff value for logFC when selecting differentially expressed genes. By default 2
#' @param pca_scale a logical variable determining whether to normalize rows when plotting PCA
#' @export
#' @return A list containing a ggplot object and a prcomp object
#' @examples
#' buildScaffold(eset_dmap,"cell_types")
#' buildScaffold(eset_dmap,"cell_types", pval_cutoff=0.01,pca_scale=T)

buildScaffold <- function(eset_scaffold,group_scaffold,pval_cutoff=0.05,lfc_cutoff=2, title="Scaffold PCA Plot",pca_scale=FALSE){
        # find DE genes
        DE_genes <- find_de_genes(eset_dmap,group_scaffold,pval_cutoff,lfc_cutoff)

        # subset
        eset_scaffold <- eset_scaffold[DE_genes,]

        # rank
        scaffold_rank<- apply(exprs(eset_scaffold),2,rank)

        # PCA
        pca <- prcomp(t(scaffold_rank),scale=pca_scale)

        # Prepare a dataframe for ggplot2
        PC1 <- pca$x[,1]
        PC2 <- pca$x[,2]
        scaffold_group <- pData(eset_scaffold)[[group_scaffold]]
        df <- data.frame(PC1,PC2,scaffold_group)

        # ggplot2
        g <- ggplot()+geom_point(data=df,mapping=aes(PC1,PC2,color=scaffold_group))+ggtitle(title)
        print(g)
        return(list(g,pca))
}
