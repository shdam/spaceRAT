#' Prefilter counts data
#'
#' @inheritParams buildScaffold
#' @noRd
preFilter <- function(object, assay = "counts", data = "logged",
                      threshold = 10) {

    # Access the counts from the SummarizedExperiment
    counts <- assay(object, assay = assay)

    # Determine genes to keep based on threshold
    if (data == "logged") {
        idx <- which(rowSums(expm1(counts)) >= threshold)
    } else if (data == "raw") {
        if (any(counts < 0)) {
            stop("Negative values are not allowed in raw count matrix!")
        }
        idx <- which(rowSums(counts) >= threshold)
    } else {
        stop("Invalid 'data' argument. Please choose 'logged' or 'raw'.")
    }

    # Check for low quality data
    if (length(idx) == 0) {
        stop("Low quality data! All genes have total counts less than ",
             threshold)
    }

    # Filter the SummarizedExperiment object
    object <- object[idx,]

    # If data is 'raw', update the assay data
    if (data == "raw") {
        assay(object, assay = assay) <- log1p(assay(object, assay = assay))
    }

    return(object)
}

