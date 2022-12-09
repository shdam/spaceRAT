#' Returns either a prebuilt scaffoldSpace object, or a new one calculated from input expression matrix
#'
#' This function returns a \code{\link{scaffoldSpace}} object that contains all parameters needed for plotting the PCA transformed sample space.
#' To use prebuilt scaffold space, simply call: buildScaffold("prebuilt_NAME"), e.g.buildScaffold("prebuilt_DMAP").
#'
#' To build your own scaffold space, pass as arguments a count matrix, a phenotype table, an a column name of the phenotype table to the function.
#' Missing values are allowed and automatically processed by the function. But please make sure that the all missing values are represented by NA, not "", " "(blank space), "-" or any other character.
#'
#' When building a user-defined scaffold space, this function first preprocesses count matrix by removing genes with total counts less than 10.
#' Then it performs differential expression analysis to select differentially expressed genes (DE genes),
#' subsets the scaffold dataset to contain only the DE genes, ranks the genes within each sample, finally performs principle component analysis (PCA) using ranks.
#'
#' By default this function plots the resulting \code{\link{scaffoldSpace}} object before returning it, but this autoplot mode can be turned off by specifying \code{auto_plot=FALSE}.
#'
#'
#' @param counts_scaffold An expression matrix of class matrix, or a data frame that can be converted to matrix. Column names are sample names.
#' Row names are gene names, which can be Ensembl Gene ID, HGNC Symbol, Entrez Gene ID Ensembl, Transcript ID or Refseq mRNA.
#' All gene ID will be automatically converted to Ensembl Gene ID. If you want to retain your gene identifier that is not currently supported, please specify \code{annotation=NA}. See parameter \code{annotation} for more information.
#' Counts of several transcript ID corresponding to same gene will be added and recorded as the counts of the gene.
#' To use the pre-built scaffold, simply set count_scaffold="prebuilt_DMAP", or "prebuilt_GTEX", and no need to specify any other parameter.

#' @param pheno_scaffold A phenotype table corresponding to the expression matrix.
#' Row names are sample names, identical to column names of \code{counts_scaffold}.
#'
#' @param colname A column name of \code{pheno_scaffold}. Cells will be grouped by this column of values.
#' Thus, differential expression analysis will be performed using this column of phenotype as independent variables.
#'
#' @param data A character indicating whether the count matrix is log-transformed or raw. By default \code{data="logged"}.
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
#' @param auto_plot A logical variable determining whether to plot the resulting scaffold space when calling the function.
#' @param annotation Type of gene identifier to use for scaffold. Currently "ensembl_gene", "ensembl_transcript", "entrez", "hgnc_symbol", and "refseq_mrna" are supported.
#' @param classes Cell types to use in plot
#' For example, set \code{annotation="hgnc_symbol"} will convert the row names (gene identifiers) of \code{counts_scaffold} to hgnc symbol, so will the \code{ExpressionSet} object and the resulting PCA scaffold.
#' If this attempted translation fails, or your desired gene identifier is not supported (especially when you are analyzing non-human data), please set \code{annotation="hgnc_symbol"} to avoid translation.
#' In this case, please manually make sure that the row names (gene identifiers) of \code{counts_scaffold} and \code{counts_sample} are the same.
#'
#' @importFrom methods is
#' @export
#' @return A scaffoldSpace object
#' @examples
#' buildScaffold("prebuilt_DMAP")
#' buildScaffold(exprs_dmap,pData_dmap,"cell_types", pval_cutoff=0.01,pca_scale=TRUE)

buildScaffold <- function(counts_scaffold,
                          pheno_scaffold = NULL,
                          colname = NULL,
                          data = "logged",
                          pcs = c(1,2),
                          plot_mode = "dot",
                          classes = NULL,
                          pval_cutoff = 0.05,
                          lfc_cutoff = 2,
                          title = "Scaffold PCA Plot",
                          pca_scale = FALSE,
                          auto_plot = TRUE,
                          annotation = "ensembl_gene"){
        # prebuilt_DMAP no samples removed
        if(is(counts_scaffold, "character") && counts_scaffold == "prebuilt_DMAP" && is(classes, "NULL")){
          space <- DMAP_scaffold
          space@pcs <- pcs
          space@plot_mode <- plot_mode
          if (auto_plot){
            g <- plotScaffold(space, title)
            print(g)
            }
          return(space)
          # prebuilt DMAP samples removed
        } else if(is(counts_scaffold, "character") && counts_scaffold == "prebuilt_DMAP" && !is(classes, "NULL")){
            se_scaffold <- createSE(counts = exprs_dmap, colData = pData_dmap, colname = "cell_types", classes = classes, annotation = annotation)
            # eset_scaffold <- createEset(exprs_dmap, pData_dmap, colname = "cell_types", classes = classes, to = annotation)
        }


        # prebuilt_GTEX missing


        # Streamline input to SummarizedExperiment
        # if(!is(counts_scaffold, "SummarizedExperiment")){
        #   eset_scaffold <- createEset(counts_scaffold, pheno_scaffold, colname, classes=classes, to=annotation)
        # }

        # data preprocessing and create eset
        if (!is(counts_scaffold, "character") && data == "logged"){
          # remove genes with total count<10
          idx <- which(rowSums(exp(counts_scaffold))<10)
          if (length(idx)==dim(counts_scaffold)[1]) stop("Low quality data! All genes have total counts less than 10.")
          if (length(idx)>0) counts_scaffold <- counts_scaffold[-idx,]
          se_scaffold <- createSE(counts_scaffold, pheno_scaffold, colname, classes = classes, annotation = annotation)
        }
        if (!is.character(counts_scaffold) && data=="raw"){
          # ensure no negative value
          if (any(counts_scaffold)<0) stop("Negative values are not allowed in raw count matrix!")
          # remove genes with total count<10
          idx <- which(rowSums(counts_scaffold)<10)
          if (length(idx)==dim(counts_scaffold)[1]) stop("Low quality data! All genes have total counts less than 10.")
          if (length(idx)>0) counts_scaffold <- counts_scaffold[-idx,]
          counts_scaffold <- log(counts_scaffold+1)
          se_scaffold <- createSE(counts_scaffold, pheno_scaffold, colname, classes = classes, annotation = annotation)

        }

        # find DE genes
        DEgenes <- findDEGenes(se_scaffold, pval_cutoff, lfc_cutoff, colname = colname)

        # subset
        se_scaffold <- se_scaffold[DEgenes,]

        # rank
        scaffold_rank <- apply(SummarizedExperiment::assay(se_scaffold, "counts"), 2, rank)

        # PCA
        pca <- stats::prcomp(t(scaffold_rank), scale = pca_scale)

        # record standard data in scaffoldSpace class
        space <- methods::new("scaffoldSpace",
                              DEgene = DEgenes,
                              label = as.character(SummarizedExperiment::colData(se_scaffold)[, colname]),
                              pca = pca,
                              pcs = pcs,
                              plot_mode = plot_mode)

        if (auto_plot){
          g <- plotScaffold(space,title)
          print(g)
        }

        return(space)
}
