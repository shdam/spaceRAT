#' Prefilter counts data
#'
#' @inheritParams buildScaffold
#' @noRd
preFilter <- function(
        mat, data,
        threshold = 10) {

    # Determine genes to keep based on threshold
    if (data == "exprs" | data == "logged") {
        idx <- which(rowSums(exp(mat)) >= threshold)
    } else if (data == "counts") {
        if (any(mat < 0)) {
            stop("Negative values are not allowed in raw count matrix!")
        }
        idx <- which(rowSums(mat) >= threshold)
    } else {
        stop("Invalid 'data' argument. Please choose 'exprs' or 'counts'.")
    }

    # Check for low quality data
    if (length(idx) == 0) {
        stop(
            "Low quality data! All genes have total counts less than ",
            threshold)
    }

    # Filter the SummarizedExperiment object
    mat <- mat[idx,]

    # If data is 'raw', convert to expression
    if (data == "counts") {
        mat <- log1p(mat)
    }

    return(mat)
}

