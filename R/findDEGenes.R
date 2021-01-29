#' Find Differentially Expressed Genes
#'
#' This function performs DE analysis using limma to find DE genes
#'
#' @importFrom Biobase pData
#' @importFrom dplyr select
#' @importFrom puma createContrastMatrix
#' @import limma
#' @param eset an ExpressionSet
#' @param group a column name of pData(eset)
#' @param pval_cutoff a cutoff value for p value when selecting differentially expressed genes. By default 0.05
#' @param lfc_cutoff a cutoff value for logFC when selecting differentially expressed genes. By default 2
#' @export
#' @return A vector of names of differentially expressed genes
#' @examples
#' find_de_genes(eset_dmap,"cell_types")

findDEGenes <- function(eset,group,pval_cutoff=0.05,lfc_cutoff=2){
  Biobase::pData(eset) <- dplyr::select(Biobase::pData(eset),group)
  design <- model.matrix(~0+Biobase::pData(eset)[[group]])
  cm <- puma::createContrastMatrix(eset) # This function only takes in ExpressionSet object
  vs_others_cols <- grep("_vs_others", colnames(cm))
  cm <- cm[,vs_others_cols]
  n <- length(unique(Biobase::pData(eset)[[group]]))
  cm[cm== -1] <- -1/(n-1)
  fit <- limma::lmFit(eset,design)
  fit <- suppressWarnings(limma::contrasts.fit(fit, contrast=cm))
  fit <- limma::eBayes(fit)
  de_genes <- sapply(colnames(cm), function(name){
    rownames(limma::topTable(fit, coef=name, number = Inf, lfc =lfc_cutoff, p.value = pval_cutoff, adjust.method="fdr", sort.by="p"))
  })
  unique(unlist(de_genes))

}
