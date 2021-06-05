#' Returns either a prebuilt scaffoldSpace object, or a new one calculated from input expression matrix
#'
#' This function returns a \code{\link{scaffoldSpace}} object that contains all parameters needed for plotting the PCA transformed sample space.
#' By default this function plots the resulting scaffoldSpace object before returning it, but this autoplot mode can be turned off by specifying \code{auto_plot=FALSE}.
#' To build your own scaffold space, pass as arguments a count matrix, a phenotype table, an a column name of the phenotype table to the function.
#' To use prebuilt scaffold space, simply call: buildScaffold("prebuilt_NAME"), e.g.buildScaffold("prebuilt_DMAP").
#'
#' When building a user-defined scaffold space, this function first preprocesses count matrix by removing genes with total counts less than 10. see \code{data} for more detail.
#' Then it performs differential expression analysis to select differentially expressed genes (DE genes),
#' subsets the scaffold dataset to contain only the DE genes, ranks the genes within each sample, finally perform principle component analysis (PCA) using ranks.
#'
#' @importFrom Biobase pData exprs
#'
#' @param counts_scaffold An expression matrix of class matrix, or a data frame that can be converted to matrix. Column names are sample names.
#' Row names are gene names, which can be Ensembl Gene ID, HGNC Symbol, Entrez Gene ID Ensembl, Transcript ID or Refseq mRNA.
#' All gene ID will be automatically converted to Ensembl Gene ID .
#' Counts of several transcript ID corresponding to same gene will be added and recorded as the counts of the gene.
#' To use the prebuilt scaffold, simply set count_scaffold="prebuilt_DMAP", or "prebuilt_GTEX", and no need to specify any other parameter.

#' @param pheno_scaffold A phenotype table corresponding to the expression matrix.
#' Row names are sample names, identical to column names of \code{counts_scaffold}.
#'
#' @param colname A column name of pheno_scaffold. Differential expression analysis will be performed using this column of phenotype as independent variables.
#' @param data A character indicating whether the count matrix is log-tansformed or raw. By default \code{data="logged"}.
#' If the count matrix is raw counts, please specify \code{data="raw"}
#'
#' @param pcs A numeric vector containing 2 numbers, indicating which two principle components to plot.
#' @param plot_mode A character indicating whether to add tiny labels to each data point.
#' By default \code{plot_mode="dot"} and tiny labels will not be attached.
#' If more than 12 cell types are to be displayed, setting \code{plot_mode="tiny_label"} may yield better visualization.
#' Shorter names for phenotypes (e.g. cell types) is strongly recommended in "tiny_label" mode.
#'
#' @param pval_cutoff A cutoff value for p value when selecting differentially expressed genes. By default \code{pval_cutoff=0.05}.
#' @param lfc_cutoff A cutoff value for logFC when selecting differentially expressed genes. By default \code{lfc_cutoff=2}.
#' @param title Title for the plot
#' @param pca_scale A logical variable determining whether to normalize rows when plotting PCA
#' @param auto_plot A logical variabe deterining whether to plot the resulting scaffold space when calling the function.
#' @export
#' @return A scaffoldSpace object
#' @examples
#' buildScaffold("prebuilt_DMAP")
#' buildScaffold("prebuilt_GTEX")
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
                          pca_scale=FALSE,
                          auto_plot=TRUE){
        # prebuilt_DMAP
        if(is.character(counts_scaffold) && counts_scaffold=="prebuilt_DMAP"){
                space<- DMAP_scaffold
                space@pcs <- pcs
                space@plot_mode <- plot_mode
                if (auto_plot){
                        g <- plotScaffold(space,title,classes=classes)
                        print(g)
                }
                return(space)
        }

        # prebuilt_GTEX missing


        # data preprocessing and create eset
        if (data=="logged"){
                # remove genes with total count<10
                idx <- which(rowSums(exp(counts_scaffold))<10)
                if (length(idx)==dim(counts_scaffold)[1]) stop("Low quality data! All genes have total counts less than 10.")
                if (length(idx)>0) counts_scaffold <- counts_scaffold[-idx,]
                eset_scaffold <- createEset(counts_scaffold,pheno_scaffold,colname)
        }
        if (data=="raw"){
                # ensure no negative value
                if (any(counts_scaffold)<0) stop("Negative values are not allowed in raw count matrix!")
                # remove genes with total count<10
                idx <- which(rowSums(counts_scaffold)<10)
                if (length(idx)==dim(counts_scaffold)[1]) stop("Low quality data! All genes have total counts less than 10.")
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
        if (auto_plot){
                g <- plotScaffold(space,title,classes=classes)
                print(g)
        }

        return(space)
}
