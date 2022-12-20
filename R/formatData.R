#' Format pheno data
#'
#' @inheritParams buildScaffold
#' @noRd
formatPheno <- function(pheno, colname = NULL, classes = NULL){

    pheno <- as.data.frame(pheno, row.names = rownames(pheno))
    # If no rownames in pheno data
    if(all( as.character(seq_len(6) %in% rownames(pheno)[seq_len(6)]))){
        pheno <- as.data.frame(pheno, row.names = pheno[[colnames(pheno)[1]]])
    }
    # If colname not in pheno
    if(is(colname, "NULL")) colname <- colnames(pheno)[1]
    if(!(colname %in% colnames(pheno)))
        stop("Column ", colname, " was not found in annotation data.")

    # Remove everything but the colname
    pheno <- pheno[, colname, drop = FALSE]

    # subset data if classes != NULL
    if (!is(classes, "NULL")){
        idx <- which(pheno[[colname]] %in% classes)
        pheno[[colname]] <- pheno[idx,,drop = FALSE]
    }

    # remove rows with NA in pheno and throw message
    complete_idx <- which(stats::complete.cases(pheno[[colname]]))
    na_rownum <- nrow(pheno)-length(complete_idx)
    if (na_rownum){
        pheno <- pheno[complete_idx,,drop = FALSE]
        message(na_rownum, " row(s) in phenotype table contain NA
                in the required column, thus removed.")
    }


    pheno[[colname]] <- as.factor(pheno[[colname]])

    return(pheno)
}

#' Remove missing annotations
#' @inheritParams buildScaffold
#' @noRd
missingAnnotation <- function(counts, pheno){
    # subset counts to only contain samples with annotation
    idx <- which(!(colnames(counts) %in% rownames(pheno)))
    if (length(idx)>0){
        counts <- counts[,-idx]
        message(length(idx),
                " sample(s) have no phenotype annotation, therefore removed.")
    }
    return(counts)

}

#' Remove annotations with missing annotations
#' @inheritParams buildScaffold
#' @noRd
missingExpression <- function(pheno, counts){
    # subset counts to only contain samples with annotation
    idx <- which(!(rownames(pheno) %in% colnames(counts)))
    if (length(idx)>0){
        pheno <- pheno[-idx,,drop = FALSE]
        message(length(idx),
                " annotation(s) have no expression data, therefore removed.")
    }
    return(pheno)
}

#' Match annotation to expression data
#' @inheritParams buildScaffold
#' @noRd
matchToExpression <- function(pheno, counts){
    # Match row/col names
    idx <- match(colnames(counts), rownames(pheno), nomatch = 0)
    idx <- idx[idx > 0]
    pheno <- pheno[idx,, drop = FALSE]
    return(pheno)
}

#' Remove NAs from expression data
#' @inheritParams buildScaffold
#' @noRd
removeNAs <- function(counts){

    # fill NA in count matrix with 0.
    is_nans <- is.na(counts)
    if (sum(is_nans) > 0){
    counts[is_nans] <- 0
    message("Expression data has ", sum(is_nans),
            " missing values. Replacing NAs by 0.")
    }
    return(counts)
}
