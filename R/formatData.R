#' Format pheno data
#'
#' @inheritParams buildScaffold
#' @importFrom stats complete.cases
#' @noRd
formatPheno <- function(object, colname = NULL, classes = NULL) {

    # Access the colData of the SummarizedExperiment object
    pheno_data <- as.data.frame(colData(object), row.names = rownames(colData(object)))

    # If no rownames in pheno_data
    if(all( as.character(seq_len(6)) %in% rownames(pheno_data)[seq_len(6)])){
        pheno_data <- as.data.frame(pheno_data, row.names = pheno_data[[colnames(pheno_data)[1]]])
    }

    # If colname not in pheno_data
    if(is.null(colname)) colname <- colnames(pheno_data)[1]
    if(!(colname %in% colnames(pheno_data)))
        stop("Column ", colname, " was not found in annotation data.")

    # Remove everything but the colname
    pheno_data <- pheno_data[, colname, drop = FALSE]

    # subset data if classes != NULL
    if (!is.null(classes)) {
        idx <- which(pheno_data[[colname]] %in% classes)
        pheno_data <- pheno_data[idx,,drop = FALSE]
    }

    # Remove rows with NA in pheno_data and throw message
    complete_idx <- which(stats::complete.cases(pheno_data[[colname]]))
    na_rownum <- nrow(pheno_data) - length(complete_idx)
    if (na_rownum) {
        pheno_data <- pheno_data[complete_idx,,drop = FALSE]
        message(na_rownum, " row(s) in phenotype table contain NA
                in the required column, thus removed.")
    }

    pheno_data[[colname]] <- factor(pheno_data[[colname]])

    # Update the colData of the SummarizedExperiment object
    object <- object[, rownames(pheno_data)]
    colData(object) <- DataFrame(pheno_data)

    return(object)
}


#' Remove missing annotations
#' @inheritParams buildScaffold
#' @noRd
missingAnnotation <- function(object) {

    # Identify samples without phenotype annotation
    idx <- which(!(colnames(object) %in% rownames(colData(object))))

    if (length(idx) > 0) {
        # Remove those samples from the SummarizedExperiment object
        object <- object[, -idx, drop = FALSE]
        message(length(idx),
                " sample(s) have no phenotype annotation, therefore removed.")
    }
    return(object)
}


#' Remove annotations with missing annotations
#' @inheritParams buildScaffold
#' @noRd
missingExpression <- function(object) {

    # Identify annotations without expression data
    idx <- which(!(rownames(colData(object)) %in% colnames(object)))

    if (length(idx) > 0) {
        # Remove those samples from the colData of the SummarizedExperiment object
        object <- object[, -idx, drop = FALSE]
        message(length(idx),
                " annotation(s) have no expression data, therefore removed.")
    }

    stopifnot(
        "None of the samples in the annotation data was found in the
        counts matrix. Please ensure the rownames and column names match"
        = ncol(object) > 0)

    return(object)
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
removeNAs <- function(object, assay = "counts") {

    # Access the counts from the SummarizedExperiment object
    counts <- assay(object, assay = assay)

    # Fill NA in count matrix with 0
    is_nans <- is.na(counts)
    if (sum(is_nans) > 0) {
        counts[is_nans] <- 0
        assay(object, assay = assay) <- counts
        message("Expression data has ", sum(is_nans),
                " missing values. Replacing NAs by 0.")
    }

    return(object)
}

