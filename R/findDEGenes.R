#' Find Differentially Expressed Genes
#'
#' This function performs DE analysis using limma to find differentially expressed genes between cell X and non-X for all types of cells that has at least 2 samples in the count matrix.
#'
#' @param eset An \code{\link{ExpressionSet}} object
#' @param pval_cutoff A cutoff value for p value when selecting differentially expressed genes. By default \code{pval_cutoff=0.05}
#' @param lfc_cutoff A cutoff value for log fold change when selecting differentially expressed genes. By default \code{lfc_cutoff=2}
#' @return A vector of names of differentially expressed genes
#' @noRd
#' @examples
#' find_de_genes(eset_dmap,"cell_types")

findDEGenes <- function(counts, cell_types, pval_cutoff = 0.05, lfc_cutoff = 2){

  message("Finding differentially expressed genes")

  # remove cell types with less than 2 cells
  count_table <- data.frame(table(cell_types))
  keep_cell_types <- count_table[count_table$Freq>1, 1]
  keep <- cell_types %in% keep_cell_types
  counts <- counts[, keep]
  cell_types <- cell_types[keep]

  # create contrast matrix
  n <- length(keep_cell_types)
  cm <- matrix(-1/(n-1),n,n)
  diag(cm) <- 1
  rownames(cm) <- keep_cell_types
  colnames(cm) <- paste(keep_cell_types,"others",sep="_")

  # limma
  design <- stats::model.matrix(~0+ cell_types)
  colnames(design) <- gsub("cell_types","",colnames(design))
  fit <- limma::lmFit(counts,design)
  fit <- limma::contrasts.fit(fit, contrast=cm)
  fit <- limma::eBayes(fit)
  de_genes <- sapply(colnames(cm), function(name){
          rownames(limma::topTable(fit, coef=name, number = Inf, lfc =lfc_cutoff, p.value = pval_cutoff, adjust.method="fdr"))
  })

  return(unique(unlist(de_genes)))

}
