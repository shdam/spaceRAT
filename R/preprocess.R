#' Preprocess data
#' @inheritParams buildScaffold
#' @noRd
preprocess <- function(
        object, colname = NULL, pheno = NULL,
        assay = "counts",
        data = "raw", threshold = 10,
        annotation = "ensembl_gene", classes = NULL){
    if(!is(object, "matrix") && !is(object, "data.frame")){
        object <- checkObject(object, assay = assay)
    } else{
        # Convert data.frame and matrix to SE
        object <- matrixToSE(object)
        if(!is.null(pheno)) colData(object) <- DataFrame(pheno)
    }
    # If colname not in colData
    if(
        !is(colname, "NULL") && !is(colData(object), "NULL") &&
        !(colname %in% names(colData(object)))) stop(
            "Column ", colname, " was not found in pheno/colData data.")
    # Convert gene names
    if (is.na(annotation) || is(annotation, "NULL")) {
        warning(
        "Gene identifier has not been resolved.
        Please manually make sure that projected samples uses same
        gene identifier as scaffold dataset")
    } else { object <- convertGeneName(object, to = annotation)}

    # Prefiltering ----
    if(is.numeric(threshold)) object <- preFilter(object, assay = assay, data, threshold)

    # Match counts and pheno ----
    if(ncol(colData(object)) > 0){
        if(is(colname, "NULL")) colname <- names(colData(object))[1]
        object <- formatPheno(object, colname, classes = classes)

        # Remove missing information
        object <- missingAnnotation(object)
        object <- missingExpression(object)

        # Ensure annotation data order matches expression data
        metadata(object)$label <- colData(object)[
                match(colnames(object), rownames(colData(object))), colname]


    } else{
        metadata(object)$label <- NULL
    }
    # Remove NAs ----
    object <- removeNAs(object, assay = assay)
    message("Preprocessing complete.")
    return(object)
}
