#' Check and extract from object
#' @inheritParams buildScaffold
#' @noRd
checkObject <- function(object, assay = "counts"){

    if(!is(object, "matrix") && !is(object, "data.frame")){
    # Convert Bioconductor objects to matrix
        if(is(object, "SummarizedExperiment")){
            pheno <- as.data.frame(object@colData)
            if(
                is(assay, "NULL") ||
                !(assay %in% names(object@assays))
                )
                assay <- names(object@assays)[1]
            object <- object@assays@data[[assay]]

        } else{stop("Expression data was not provided in a supported format
                    (matrix, data.frame, SummarizedExperiment,
                    or SingleCellExperiment).\n If you'd like us to
                    support a new format, please raise an issue on GitHub.")}
        return(list(object, pheno))
        } else{
            return(object)
        }
    }
