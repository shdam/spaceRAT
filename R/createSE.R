#' Create SummarizedExperiment object
#'
#' This function creates SummarizedExperiment from an expression matrix and a phenotype data frame.
#' This function also performs data preprocessing: samples lacking either count data or phenotype annotation are removed. NAs are replaced by 0.
#'
#' @param counts An expression matrix of class matrix or data frame that can be converted to matrix. Each row corresponds to a gene. Each column corresponds to a sample.
#' @param colData A dataframe of sample annotation corresponding to the expression matrix. Each row corresponds to a sample. Each column corresponds to a sample feature.
#' @param colname A column name of \code{pheno}. Differential expression analysis will be performed using this column of phenotype as independent variables.
#' @return A SummarizedExperiment object
#' @noRd
#' @examples
#' createSE(exprs_dmap[1:20,1:10],pData_dmap[1:10,,drop=FALSE])
createSE <- function(counts, colData, colname = "cancer_type", classes = NULL, to = "ensembl_gene"){

  # If counts is a data.frame, check for non-numeric values
  if(is(counts, "data.frame")){
    numerics <- sapply(counts, is, "numeric")
    stopifnot("To many non-numerics in counts data! Please separate from counts." = sum(!numerics) == 1)
    # Only one non-numeric - use as rownames
    counts <- tibble::column_to_rownames(counts, colnames(counts)[which(!numerics)]) %>%
      as.matrix()
  }

  # check class of counts, convert to matrix
  if (!is(counts, "matrix")){
    tryCatch(counts <- data.matrix(counts),
             error = function(e) {
               stop("Expression matrix cannot be converted to matrix:", e)
             })
  }

  # check class of colData, convert to dataframe
  if (!is(colData, "NULL")){
    # If colname not in colData
    if(!(colname %in% colnames(colData))) stop(paste("Column", colname, "was not found in annotation data."))
    tryCatch(colData <- S4Vectors::DataFrame(colData, row.names = colData[[colname]]),
             error=function(e){
               stop("Phenotype data cannot be converted to DataFrame:", e)
             })
    colnames(counts) <- rownames(colData)
  }

  # convert gene names
  if (is.na(to) || is(to, "NULL")) {
    warning("Gene identifier has not been resolved. Please manually make sure that projected samples uses same gene identifier as scaffold dataset")
  } else {
    counts <- convertGeneName(counts, to = to)
    }

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list("counts" = counts),
    colData = colData
  )

  # remove rows with NA in colData and throw message
  complete_idx <- which(stats::complete.cases(colnames(se)))
  na_rownum <- ncol(se)-length(complete_idx)
  if (na_rownum){
    se <- se[, complete_idx, drop = FALSE]
    message(na_rownum, " row(s) in phenotype table contain NA in the required column, it is therefore removed.")
  }

  # subset counts to only contain samples with annotation
  idx <- which(!(colnames(counts) %in% rownames(pheno)))
  if (length(idx)>0){
    counts <- counts[,-idx,drop=F]
    message(length(idx)," sample(s) have no phenotype annotation, thus removed.")
  }
  idx <- which(! rownames(pheno) %in% colnames(counts))
  if (length(idx)>0){
    pheno <- pheno[,-idx,drop=F]
    message(length(idx)," sample(s) have no expression count data, thus removed.")
  }

  # match row/col names
  idx <- match(colnames(counts),rownames(pheno))
  pheno <- pheno[idx,,drop=F]

  # fill NA in count matrix with 0.
  sum_na <- sum(is.na(counts))
  if (sum_na>0){
    counts[is.na(counts)] <- 0
    message("Count matrix has ",sum_na, " missing values. Replace NA by 0.")
  }

  # create ExpressionSet
  eset <- Biobase::ExpressionSet(counts, phenoData = Biobase::AnnotatedDataFrame(pheno))

  # subset ExpressionSet if classes!=NULL
  if (!is.null(classes)){
    idx <- which(Biobase::pData(eset)[[colname]] %in% classes)
    eset <- eset[,idx]
    Biobase::pData(eset)[,1] <-droplevels(Biobase::pData(eset)[,1])
  }

  return(eset)



}
