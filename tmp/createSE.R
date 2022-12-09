#' Create SummarizedExperiment object
#'
#' This function creates SummarizedExperiment from an expression matrix and a phenotype data frame.
#' This function also performs data preprocessing: samples lacking either count data or phenotype annotation are removed. NAs are replaced by 0.
#'
#' @inheritParams buildScaffold
#' @param counts An expression matrix of class matrix or data frame that can be converted to matrix. Each row corresponds to a gene. Each column corresponds to a sample.
#' @param colData A dataframe of sample annotation corresponding to the expression matrix. Each row corresponds to a sample. Each column corresponds to a sample feature.
#' @param colname A column name of \code{pheno}. Differential expression analysis will be performed using this column of phenotype as independent variables.
#' @return A SummarizedExperiment object
#' @noRd
#' @examples
#' createSE(exprs_dmap[1:20,1:10],pData_dmap[1:10,,drop=FALSE])
createSE <- function(counts, colData = NULL, colname = "cancer_type", classes = NULL, annotation = "ensembl_gene"){

  ### Error handling ----

  stopifnot("Please provide expression matrix only as a matrix or data.frame" = (is(counts, "data.frame") || is(counts, "matrix")))

  # Identify non-numerics in counts data
  numerics <- apply(counts, 2, is, "numeric")
  stopifnot("Too many non-numerics in counts data! Please separate from counts." = sum(!numerics) < 2)
  empty_rownames <- all(rownames(counts)[1:6] == paste(1:6))

  # No rowname information
  stopifnot("The expression matrix does not contain rownames" = (sum(numerics) == 1 || !empty_rownames))

  # Check counts rownames
  if(sum(numerics) == 1){

    if(empty_rownames){ # No rownames
      if(is(counts, "data.frame")){
        counts <- tibble::column_to_rownames(counts, colnames(counts)[which(!numerics)]) %>%
          as.matrix()
      } else if(is(counts, "matrix")){
        rownames(counts) <- counts[which(!numerics)]
        counts <- counts[which(numerics)]
      }
    } else if(!empty_rownames){ # Rownames and additional column
      warning(paste("column:", colnames(counts)[which(!numerics)], "is removed from expression matrix. If that column was intended to be the matrix' rownames, please manually set the rownames."))
      counts <- counts[which(numerics)]
    }
  } else if(is(counts, "data.frame")) {
    counts <- as.matrix(counts)
  } # Else the matrix is fine

  ### convert gene names ----
  if (is.na(annotation) || is(annotation, "NULL")) {
    warning("Gene identifier has not been resolved. Please manually make sure that projected samples uses same gene identifier as scaffold dataset")
  } else {
    counts <- convertGeneName(counts, to = annotation)
  }


  if(is(colData, "NULL")){
    warning("No phenotype information was provided. Please make sure the sample colnames match the scaffold.")
    se <- SummarizedExperiment::SummarizedExperiment(
      assays = list("counts" = counts)
    )
  } else{

    if(!is(colData, "DFrame")){
      # If no rownames in colData
      if(all(rownames(colData)[1:6] == paste(1:6))){
        colData <- S4Vectors::DataFrame(colData, row.names = colData[[colnames(colData)[1]]])
      } else{# if rownames already present
        colData <- S4Vectors::DataFrame(colData, row.names = rownames(colData))
      }
    }
    # If colname not in colData
    if(!is(colname, "NULL") && !(colname %in% colnames(colData))) stop(paste("Column", colname, "was not found in annotation data."))

    # Ensure data matches
    idx <- match(colnames(counts), rownames(colData))
    counts <- counts[, idx]
    colData <- colData[idx,, drop = FALSE]

    se <- SummarizedExperiment::SummarizedExperiment(
      assays = list("counts" = counts),
      # assays = counts,
      colData = colData
    )
    # Should it be like this?
    # colnames(counts) <- rownames(colData)
  }






  # subset SummarizedExperiment if classes!=NULL
  if (!is.null(classes)){
    idx <- which(colnames(se) %in% classes)
    se <- se[,idx]
    # Biobase::pData(eset)[,1] <- droplevels(Biobase::pData(eset)[,1])
  }

  return(se)


  # remove rows with NA in colData and throw message
  complete_idx <- which(stats::complete.cases(colnames(se)))
  na_rownum <- ncol(se)-length(complete_idx)
  if (na_rownum > 0){
    se <- se[, complete_idx, drop = FALSE]
    message(na_rownum, " row(s) in phenotype table contain NA in the required column, it is therefore removed.")
  }

  # subset counts to only contain samples with annotation
  idx <- which(!(colnames(counts) %in% rownames(pheno)))
  if (length(idx)>0){
    counts <- counts[,-idx,drop=F]
    message(length(idx)," sample(s) have no phenotype annotation, thus removed.")
  }
  idx <- which(!(rownames(pheno) %in% colnames(counts)))
  if (length(idx)>0){
    pheno <- pheno[,-idx,drop=F]
    message(length(idx)," sample(s) have no expression count data, thus removed.")
  }

  # match row/col names
  idx <- match(colnames(counts), rownames(pheno))
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
