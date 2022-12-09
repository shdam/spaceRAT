#' Find Differentially Expressed Genes
#'
#' This function performs DE analysis using limma to find differentially expressed genes between cell X and non-X for all types of cells that has at least 2 samples in the count matrix.
#'
#' @param se A \code{\link{SummarizedExperiment}} object
#' @param pval_cutoff A cutoff value for p value when selecting differentially expressed genes. By default \code{pval_cutoff=0.05}
#' @param lfc_cutoff A cutoff value for log fold change when selecting differentially expressed genes. By default \code{lfc_cutoff=2}
#' @return A vector of names of differential expressed genes
#' @noRd
#' @examples
#' # findDEGenes(se_dmap,"cell_types")

findDEGenes_SE <- function(se, pval_cutoff = 0.05, lfc_cutoff = 2, colname = NULL){
        # remove cell types with less than 2 cells
        count_table <- data.frame(table(SummarizedExperiment::colData(se)[, colname]))
        cell_types <- count_table[count_table$Freq>1,]$Var1
        se <- se[, SummarizedExperiment::colData(se)[, colname] %in% cell_types]

        # create contrast matrix
        n <- length(cell_types)
        cm <- matrix(-1/(n-1),n,n)
        diag(cm) <- 1
        rownames(cm) <- cell_types
        colnames(cm) <- paste(cell_types, "others", sep = "_")

        # limma
        cell_types <- SummarizedExperiment::colData(se)[, colname]
        design <- stats::model.matrix(~0+ cell_types)
        colnames(design) <- gsub("cell_types","",colnames(design))
        fit <- limma::lmFit(eset,design)
        fit <- limma::contrasts.fit(fit, contrast=cm)
        fit <- limma::eBayes(fit)
        de_genes <- sapply(colnames(cm), function(name){
                rownames(limma::topTable(fit, coef=name, number = Inf, lfc =lfc_cutoff, p.value = pval_cutoff, adjust.method="fdr"))
        })

        return(unique(unlist(de_genes)))

}
