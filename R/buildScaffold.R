#' Returns a ggplot and PCA parameters for ExpressionSet object
#'
#' This function performs differential expression analysis to select DE gene, then subset the scaffold dataset, rank the subset, finally produce a PCA plot using ranks.
#' @import dplyr
#' @import magrittr
#' @import  ggplot2
#' @importFrom Biobase pData exprs
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @param counts_scaffold an expression matrix of class matrix, or data frame that can be converted to matrix. To use our prebuilt scaffold, you can set count_scaffold="prebuilt_DMAP", or "prebuilt_GTEX".
#' @param pheno_scaffold a phenotype table corresponding to the expression matrix.
#' @param colname a column name of pData_scaffold. Differential analysis will be performed using this column of phenotype as independent variables.
#' @param data a character indicating whether the count matrix is log-tansformed or raw. To perform DE analysis, dataset need to be log-transformed, so setting "data=raw" will automatically perform log transformation before DE analysis.
#' @param pcs a numeric vector containing 2 numbers, indicating which two principle components to plot.
#' @param plot_mode a character indicating whether to add tiny labels to each data point.
#' @param title title for the plot
#' @param pval_cutoff a cutoff value for p value when selecting differentially expressed genes. By default 0.05
#' @param lfc_cutoff a cutoff value for logFC when selecting differentially expressed genes. By default 2
#' @param pca_scale a logical variable determining whether to normalize rows when plotting PCA
#' @export
#' @return A scaffoldSpace object
#' @examples
#' buildScaffold(exprs_dmap,pData_dmap"cell_types")
#' buildScaffold(exprs_dmap,pData_dmap,"cell_types", pval_cutoff=0.01,pca_scale=T)

buildScaffold <- function(counts_scaffold,
                          pheno_scaffold,
                          colname,
                          data="logged",
                          pcs=c(1,2),
                          plot_mode="dot",
                          classes = NULL,
                          pval_cutoff=0.05,
                          lfc_cutoff=2,
                          title="Scaffold PCA Plot",
                          pca_scale=FALSE){

        if(is.character(counts_scaffold) && counts_scaffold=="prebuilt_DMAP"){
                space<- DMAP_scaffold
                space@pcs <- pcs
                space@plot_mode <- plot_mode
                g <- plotScaffold(space,title,classes=classes)
                print(g)
                return(space)
        }

        # data preprocessing and create eset
        if (data=="logged"){
                # remove genes with total count<10
                idx <- which(rowSums(exp(counts_scaffold))<10)
                if (length(idx)>0) counts_scaffold <- counts_scaffold[-idx,]
                eset_scaffold <- createEset(counts_scaffold,pheno_scaffold,colname)
        }
        if (data=="raw"){
                # remove genes with total count<10
                idx <- which(rowSums(counts_scaffold)<10)
                if (length(idx)>0) counts_scaffold <- counts_scaffold[-idx,]
                counts_scaffold <- log(counts_scafold+1)
                eset_scaffold <- createEset(counts_scaffold,pheno_scaffold,colname)
        }

        # find DE genes
        DEgenes <- findDEGenes(eset_scaffold,pval_cutoff,lfc_cutoff)

        # subset
        eset_scaffold <- eset_scaffold[DEgenes,]

        # rank
        scaffold_rank<- apply(Biobase::exprs(eset_scaffold),2,rank)

        # PCA
        pca <- prcomp(t(scaffold_rank),scale=pca_scale)

        # record standard data in scaffoldSpace class
        space <- new("scaffoldSpace",
                     DEgene=DEgenes,
                     label=as.character(Biobase::pData(eset_scaffold)[[colname]]),
                     pca=pca,
                     pcs=pcs,
                     plot_mode=plot_mode)

        g <- plotScaffold(space,title,classes=classes)
        print(g)

        return(space)
}
