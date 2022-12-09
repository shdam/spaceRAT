preprocess <- function(object,
                       colname = NULL,
                       pheno = NULL,
                       assay = "counts",
                       data = "logged",
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


  # Convert gene names
  if (is.na(annotation) || is(annotation, "NULL")) {
    warning("Gene identifier has not been resolved. Please manually make sure that projected samples uses same gene identifier as scaffold dataset")
  } else { counts <- convertGeneName(counts, to = annotation)}

  # (missing) Prefiltering ----
  # data preprocessing and create eset
  if (data == "logged"){
    # remove genes with total count<10
    idx <- which(rowSums(exp(counts))<10)
    if (length(idx)==dim(counts)[1]) stop("Low quality data! All genes have total counts less than 10.")
    if (length(idx)>0) counts <- counts[-idx,]

  } else if(data=="raw"){
    # ensure no negative value
    if (any(counts<0)) stop("Negative values are not allowed in raw count matrix!")
    # remove genes with total count<10
    idx <- which(rowSums(counts)<10)
    if (length(idx)==dim(counts)[1]) stop("Low quality data! All genes have total counts less than 10.")
    if (length(idx)>0) counts <- counts[-idx,]
    counts <- log(counts+1)
  }

  # Match counts and pheno ----
  if(!is(pheno, "NULL")){
    if(is(colname, "NULL")) colname <- colnames(pheno)[1]
    pheno <- formatPheno(pheno, colname)

    # Remove missing information
    counts <- missingAnnotation(counts, pheno)
    pheno <- missingExpression(pheno, counts)

    # Ensure annotation data order matches expression data
    cell_types <- pheno[match(colnames(counts), rownames(pheno)), ]

  }


  # Remove NAs
  counts <- removeNAs(counts)
  message("Preprocessing complete.")

  return(list(counts, cell_types))
}
