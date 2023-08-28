#' Format pheno data
#'
#' @inheritParams buildScaffold
#' @importFrom stats complete.cases
#' @noRd
formatPheno <- function(pheno, colname = NULL, classes = NULL) {

    # If no rownames in pheno_data
    if(all( as.character(seq_len(6)) %in% rownames(pheno)[seq_len(6)])){
        pheno <- as.data.frame(pheno[, -1, drop=FALSE], row.names = as.character(pheno[[colnames(pheno)[1]]]))
    }

    # If colname not in pheno_data
    if(is.null(colname)) colname <- colnames(pheno)[1]
    if(!(colname %in% colnames(pheno)))
        stop("Column ", colname, " was not found in annotation data.")

    # Remove everything but the colname
    pheno <- pheno[, colname, drop = FALSE]

    # subset data if classes != NULL
    if (!is.null(classes)) {
        idx <- which(pheno[[colname]] %in% classes)
        pheno <- pheno[idx,,drop = FALSE]
    }

    # Remove rows with NA in pheno_data and throw message
    complete_idx <- which(stats::complete.cases(pheno[[colname]]))
    na_rownum <- nrow(pheno) - length(complete_idx)
    if (na_rownum) {
        pheno <- pheno[complete_idx,,drop = FALSE]
        message(na_rownum, " row(s) in phenotype table contain NA
                in the required column, thus removed.")
    }

    pheno[[colname]] <- factor(pheno[[colname]])
    return(pheno)
}


#' Remove missing annotations
#' @inheritParams buildScaffold
#' @noRd
missingAnnotation <- function(mat, pheno) {

    # Identify samples without phenotype annotation
    idx <- which(!(colnames(mat) %in% rownames(pheno)))

    if (length(idx) > 0) {
        # Remove those samples from the SummarizedExperiment object
        mat <- mat[, -idx, drop = FALSE]
        message(
            length(idx),
            " sample(s) have no phenotype annotation, therefore removed.")
    }
    return(mat)
}


#' Remove annotations with missing annotations
#' @inheritParams buildScaffold
#' @noRd
missingExpression <- function(mat, pheno) {

    # Identify annotations without expression data
    idx <- which(!(rownames(pheno) %in% colnames(mat)))

    if (length(idx) > 0) {
        # Remove those samples from the pheno data
        pheno <- pheno[-idx, , drop = FALSE]
        message(
            length(idx),
            " annotation(s) have no expression data, therefore removed.")
    }

    stopifnot(
        "None of the samples in the annotation data was found in the
        counts matrix. Please ensure the rownames and column names match"
        = nrow(pheno) > 0)

    return(pheno)
}


#' Match annotation to expression data
#' @inheritParams buildScaffold
#' @noRd
matchToExpression <- function(mat, pheno){
    # Match row/col names
    idx <- match(colnames(mat), rownames(pheno), nomatch = 0)
    idx <- idx[idx > 0]
    pheno <- pheno[idx,, drop = FALSE]
    return(pheno)
}

#' Remove NAs from expression data
#' @inheritParams buildScaffold
#' @noRd
removeNAs <- function(mat) {
    # Fill NA in count matrix with 0
    is_nans <- is.na(mat)
    if (sum(is_nans) > 0) {
        mat[is_nans] <- 0
        message(
            "Expression data has ", sum(is_nans),
            " missing values. Replacing NAs by 0.")
    }
    return(mat)
}

