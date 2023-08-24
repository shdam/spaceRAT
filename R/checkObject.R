#' Check and extract from object
#' @inheritParams buildScaffold
#' @noRd
checkObject <- function(object, assay = "counts") {

    if (is(object, "SummarizedExperiment")) {
        # Check if the desired assay is present
        if (assay %in% SummarizedExperiment::assayNames(object)) {
            return(object)
        } else {
            stop("The provided SummarizedExperiment does not contain a '",
                 assay, "' assay. Available assays: ",
                 paste(SummarizedExperiment::assayNames(object), collapse=", "))
        }

    } else if (is(object, "matrix") || is(object, "data.frame")) {
        return(object)
    } else {
        stop("Expression data was not provided in a supported format ",
             "(matrix, data.frame, or SummarizedExperiment).\n",
             "If you'd like us to support a new format, ",
             "please raise an issue on GitHub.")
    }
}



