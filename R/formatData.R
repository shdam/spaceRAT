#' Format pheno data
#'
#'
formatPheno <- function(pheno, colname = NULL, classes = NULL){

  pheno <- as.data.frame(pheno, row.names = rownames(pheno))
  # If no rownames in pheno data
  if(all(rownames(pheno)[1:6] == paste(1:6))){
    pheno <- as.data.frame(pheno, row.names = pheno[[colnames(pheno)[1]]])
  }
  # If colname not in pheno
  if(is(colname, "NULL")) colname <- colnames(pheno_scaffold)[1]
  if(!(colname %in% colnames(pheno))) stop(paste("Column", colname, "was not found in annotation data."))

  # Remove everything but the colname
  pheno <- pheno[, colname, drop = FALSE]

  # subset data if classes != NULL
  if (!is(classes, "NULL")){
    idx <- which(pheno[[colname]] %in% classes)
    pheno[[colname]] <- pheno[idx,,drop=F]
  }

  # remove rows with NA in pheno and throw message
  complete_idx <- which(stats::complete.cases(pheno[[colname]]))
  na_rownum <- nrow(pheno)-length(complete_idx)
  if (na_rownum){
    pheno <- pheno[complete_idx,,drop=F]
    message(na_rownum, " row(s) in phenotype table contain NA in the required column, thus removed.")
  }


  pheno[[colname]] <- as.factor(pheno[[colname]])

  return(pheno)
}


missingAnnotation <- function(counts, pheno){
  # subset counts to only contain samples with annotation
  idx <- which(!(colnames(counts) %in% rownames(pheno)))
  if (length(idx)>0){
    counts <- counts[,-idx]
    message(length(idx)," sample(s) have no phenotype annotation, therefore removed.")
  }
  return(counts)

}

missingExpression <- function(pheno, counts){
  # subset counts to only contain samples with annotation
  idx <- which(!(rownames(pheno) %in% colnames(counts)))
  if (length(idx)>0){
    pheno <- pheno[-idx,,drop=F]
    message(length(idx)," annotation(s) have no expression data, therefore removed.")
  }
  return(pheno)
}


matchToExpression <- function(pheno, counts){
  # Match row/col names
  idx <- match(colnames(counts), rownames(pheno), nomatch = 0)
  idx <- idx[idx > 0]
  pheno <- pheno[idx,, drop = FALSE]
  return(pheno)
}


removeNAs <- function(counts){

  # fill NA in count matrix with 0.
  is_nans <- is.na(counts)
  if (sum(is_nans) > 0){
    counts[is_nans] <- 0
    message("Expression data has ", sum(is_nans), " missing values. Replacing NAs by 0.")
  }
  return(counts)

}
