#' Preprocess data
#' @inheritParams buildScaffold
#' @importFrom SummarizedExperiment assay colData
#' @noRd
preprocess <- function(
        object, colname = NULL, pheno = NULL,
        assay = NULL,
        data = NULL, threshold = 10,
        annotation = "ensembl_gene", classes = NULL){
    if(!is(object, "matrix") && !is(object, "data.frame")){
        object <- checkObject(object, assay = assay)
        mat <- assay(object, assay)
        if(is(pheno, "NULL")) pheno <- as.data.frame(
            colData(object), row.names = rownames(colData(object)))
        rm(object)
    } else{# Check matrix or data.frame
        mat <- checkMatrix(object)
    }
    if( # If colname not in pheno data
        !is(colname, "NULL") && !is(pheno, "NULL") &&
        !(colname %in% colnames(pheno))) stop(
            "Column ", colname, " was not found in pheno/colData data.")
    # Convert gene names
    if (is.na(annotation) || is(annotation, "NULL")) {
        warning(
        "Gene identifier has not been resolved.
        Please manually make sure that projected samples uses same
        gene identifier as scaffold dataset")
    } else { mat <- convertGeneName(mat, to = annotation)}

    # Prefiltering ----
    if(is.numeric(threshold)) mat <- preFilter(
        mat, data = data, threshold = threshold)
    # Match counts and pheno ----
    if(!is(pheno, "NULL")){
        if(is(colname, "NULL")) colname <- colnames(pheno)[1]
        pheno <- formatPheno(pheno, colname, classes = classes)

        # Remove missing information
        mat <- missingAnnotation(mat, pheno)
        pheno <- missingExpression(mat, pheno)

        # Ensure annotation data order matches expression data
        pheno <- matchToExpression(mat, pheno)

    }
    # Remove NAs ----
    mat <- removeNAs(mat)
    message("Preprocessing complete.")
    return(list("mat" = mat, "pheno" = pheno))
}
