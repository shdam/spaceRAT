#' Convert matrix to SummarizedExperiment
#'
#'
#' @param mat A matrix to check
#' @noRd
matrixToSE <- function(mat){

    # Identify non-numerics in counts data
    numerics <- vapply(mat[seq_len(3),], is, "numeric", FUN.VALUE = TRUE)
    stopifnot(
        "Too many non-numerics in counts data! Please separate from counts." =
            sum(!numerics) < 2)
    empty_rownames <- all(
        as.character(seq_len(3)) %in% rownames(mat)[seq_len(3)]
    ) || is.null(rownames(mat))

    # No rowname information
    stopifnot(
        "The expression matrix does not contain rownames" =
            (sum(!numerics) == 1 || !empty_rownames))

    # Check counts rownames
    if(sum(!numerics) == 1){
        if(empty_rownames){ # No rownames
            if(is(mat, "data.frame")){
                mat <- as.data.frame(mat)
                rownames(mat) <- mat[[colnames(mat)[which(!numerics)]]]
                mat[[colnames(mat)[which(!numerics)]]] <- NULL
            }
        } # Else the matrix is fine
    }

    # Convert to matrix
    mat <- as.matrix(mat)

    # Create SummarizedExperiment object
    se <- SummarizedExperiment(assays = list(counts = mat))

    return(se)
}

