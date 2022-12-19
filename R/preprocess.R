#' Preprocess data
#' @inheritParams buildScaffold
#' @noRd
preprocess <- function(object,
                       colname = NULL,
                       pheno = NULL,
                       assay = "counts",
                       data = "logged",
                       threshold = 10,
                       annotation = "ensembl_gene",
                       classes = NULL){


    if(!is(object, "matrix") && !is(object, "data.frame")){
        object <- checkObject(object, assay = assay)
        counts <- object[[1]]
        pheno <- object[[2]]
    } else{
        # Check that matrix has rownames (converts data.frame to matrix)
        counts <- checkMatrix(object)
    }
    rm(object)

    # If colname not in pheno
    if(!is(colname, "NULL") && !is(pheno, "NULL") &&
       !(colname %in% colnames(pheno))) stop(
           "Column ", colname, " wasn't found in pheno data."
           )

    # Convert gene names
    if (is.na(annotation) || is(annotation, "NULL")) {
        warning("Gene identifier has not been resolved.
                Please manually make sure that projected samples uses same
                gene identifier as scaffold dataset")
    } else { counts <- convertGeneName(counts, to = annotation)}

    # Prefiltering ----
    if(is.numeric(threshold)) counts <- pre_filter(counts, data, threshold)

    # Match counts and pheno ----
    if(!is(pheno, "NULL")){
        if(is(colname, "NULL")) colname <- colnames(pheno)[1]
        pheno <- formatPheno(pheno, colname)

        # Remove missing information
        counts <- missingAnnotation(counts, pheno)
        pheno <- missingExpression(pheno, counts)

        # Ensure annotation data order matches expression data
        cell_types <- pheno[match(colnames(counts), rownames(pheno)), ]

    } else{
        cell_types <- NULL
    }


    # Remove NAs ----
    counts <- removeNAs(counts)
    message("Preprocessing complete.")

    return(list(counts, cell_types))
}
