#' Find Differentially Expressed Genes
#'
#' This function performs DE analysis using limma to find
#'   differentially expressed genes between cell X and non-X for all
#'   types of cells that has at least 2 samples in the count matrix.
#'
#' @param pval_cutoff A cutoff value for p value when selecting
#' differentially expressed genes. By default \code{pval_cutoff=0.05}
#' @param lfc_cutoff A cutoff value for log fold change when selecting
#' differentially expressed genes. By default \code{lfc_cutoff=2}
#' @return A vector of names of differentially expressed genes
#' @importFrom stats model.matrix
#' @importFrom limma lmFit contrasts.fit eBayes topTable
#' @noRd
findDEGenes <- function(mat, label, pval_cutoff = 0.05,
                        lfc_cutoff = 2){

    # remove cell types with less than 2 cells
    count_table <- data.frame(table(label))
    keep_labels <- count_table[count_table$Freq>1, 1]
    stopifnot(
        "All cells have unique phenotype information.
        Please group cells by phenotype." = length(keep_labels) > 0)
    keep <- label %in% keep_labels
    mat <- mat[, keep]
    label <- label[keep]

    # create contrast matrix
    n <- length(keep_labels)
    cm <- matrix(-1/(n-1),n,n)
    diag(cm) <- 1
    rownames(cm) <- keep_labels
    colnames(cm) <- paste(keep_labels,"others",sep="_")

    # limma
    design <- stats::model.matrix(~0+ label)
    colnames(design) <- gsub("label","", colnames(design))
    fit <- limma::lmFit(mat, design)
    fit <- limma::contrasts.fit(fit, contrast=cm)
    fit <- limma::eBayes(fit)
    de_genes <- lapply(colnames(cm), function(name){
        rownames(limma::topTable(
            fit, coef = name, number = Inf,
            lfc = lfc_cutoff, p.value = pval_cutoff,
            adjust.method = "fdr"))
    })
    names(de_genes) <- keep_labels
    return(de_genes)
}


findDEGenes2 <- function(mat, group, labels, pval_cutoff = 0.05, lfc_cutoff = 2,
                         n_genes = 15) {

    # Validate that group is in labels
    if (!(group %in% labels)) {
        stop("Specified group not found in labels.")
    }

    # Create binary vector indicating whether each label is the group of interest
    group_indicator <- as.factor(labels == group)

    # Limma analysis
    design <- model.matrix(~ group_indicator)
    fit <- limma::lmFit(mat, design)

    # Define contrast for group vs. all others
    # The contrast matrix should have two rows to match the coefficients in the linear model
    contrast_matrix <- matrix(c(0, 1), nrow = 2)  # Compare the group of interest against others
    colnames(contrast_matrix) <- "GroupEffect"

    fit <- limma::contrasts.fit(fit, contrasts = contrast_matrix)
    fit <- limma::eBayes(fit)

    # Extract differentially expressed genes for the specified contrast
    de_genes <- rownames(limma::topTable(
        fit, coef = "GroupEffect", number = n_genes,
        lfc = lfc_cutoff, p.value = pval_cutoff, sort.by = "logFC",
        adjust.method = "fdr"
    ))

    return(de_genes)
}
