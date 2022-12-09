#' Create ExpressionSet object
#'
#' This function creates ExpressionSet from an expression matrix and a phenotype data frame.
#' This function also performs data preprocessing: samples lacking either count data or phenotype annotation are removed. NAs are replaced by 0.
#'
#' @inheritParams buildScaffold
#' @param counts An expression matrix of class matrix, or data frame that can be converted to matrix. Each row corresponds to a gene. Each column corresponds to a sample.
#' @param pheno A dataframe of sample annotation corresponding to the expression matrix. Each row corresponds to a sample. Each column corresponds to a sample feature.
#' @param colname A column name of \code{pheno}. Differential expression analysis will be performed using this column of phenotype as independent variables.
#' @return An ExpressionSet object
#' @noRd
#' @examples
#' createEset(exprs_dmap[1:20,1:10],pData_dmap[1:10,,drop=FALSE])


createEset <- function(counts, pheno = NULL, colname = "cancer_type", classes = NULL, annotation = "ensembl_gene"){

  ### Error handling ----

  stopifnot("Please provide expression matrix only as a matrix or data.frame" = (is(counts, "data.frame") || is(counts, "matrix")))

  # Check counts matrix
  counts <- checkMatrix(counts)

  # Create minimal ExpressionSet
  eset <- Biobase::ExpressionSet(assayData = counts,
                                 annotation = annotation)
  # Remove counts matrix to save memory
  rm(counts)

  # Check pheno data
  if(!is(pheno, "NULL")){

    pheno <- BiocGenerics::as.data.frame(pheno, row.names = rownames(pheno))
    # If no rownames in pheno data
    if(all(rownames(pheno)[1:6] == paste(1:6))){
      pheno <- BiocGenerics::as.data.frame(pheno, row.names = pheno[[colnames(pheno)[1]]])
    }
    # If colname not in pheno
    if(!is(colname, "NULL") && !(colname %in% colnames(pheno))) stop(paste("Column", colname, "was not found in annotation data."))


  }

  # convert gene names
  if (is.na(annotation) || is(annotation, "NULL")) {
          warning("Gene identifier has not been resolved. please manually make sure that projected samples uses same gene identifier as scaffold dataset")
  } else { eset <- convertGeneName(eset, to = annotation)}



  #select the specified column of phenotype table as final phenotype table
  # pheno <- pheno[,colname,drop=F]
  # pheno[[colname]] <- as.factor(pheno[[colname]])
  # pheno[,1] <- as.factor(pheno[,1])

  # remove rows with NA in pheno and throw message
  complete_idx <- which(stats::complete.cases(pheno[[colname]]))
  na_rownum <- dim(pheno)[1]-length(complete_idx)
  if (na_rownum){
    pheno <- pheno[complete_idx,,drop=F]
    message(na_rownum, " row(s) in phenotype table contain NA in the required column, thus removed.")
  }

  # subset counts to only contain samples with annotation
  idx <- which(!(colnames(eset) %in% rownames(pheno)))
  if (length(idx)>0){
    eset <- eset[,-idx]
    message(length(idx)," sample(s) have no phenotype annotation, thus removed.")

  }
  idx <- which(!(rownames(pheno) %in% colnames(eset)))
  if (length(idx)>0){
    pheno <- pheno[-idx,,drop=F]
    message(length(idx)," sample(s) have no expression count data, thus removed.")
  }


  # match row/col names
  # Ensure data matches
  idx <- match(colnames(eset), rownames(pheno), nomatch = 0)
  idx <- idx[idx > 0]
  pheno <- pheno[idx,, drop = FALSE]

  pheno[[colname]] <- as.factor(pheno[[colname]])

  # fill NA in count matrix with 0.
  is_nans <- is.na(Biobase::exprs(eset))
  if (sum(is_nans) > 0){
    Biobase::exprs(eset)[is_nans] <- 0
    message("Count matrix has ", sum(is_nans), " missing values. Replace NA by 0.")
  }

  # create ExpressionSet
  Biobase::pData(eset) <- Biobase::pData(Biobase::AnnotatedDataFrame(pheno))

  # subset ExpressionSet if classes!=NULL
  if (!is.null(classes)){
    idx <- which(Biobase::pData(eset)[, colname] %in% classes)
    eset <- eset[,idx]
    Biobase::pData(eset)[, colname] <- droplevels(Biobase::pData(eset)[, colname])
  }

  return(eset)



}
